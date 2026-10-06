/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include <algorithm>
#include <vector>
#include "move_surf.h"
#include "sparta.h"
#include "surf.h"
#include "surf_collide.h"
#include "grid.h"
#include "particle.h"
#include "comm.h"
#include "update.h"
#include "domain.h"
#include "compute.h"
#include "geometry.h"
#include "input.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathExtra;
using namespace MathConst;

enum{READFILE,TRANSLATE,ROTATE};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

// [ADDED] Phase 1: engulf キーワードの値

enum{ENGULF_DELETE,ENGULF_PROJECT};

#define MAXLINE 256

// [ADDED] Phase 1: 退避処理の定数とデバッグフラグ
// EPSPUSH = retreat distance from surf, as fraction of local cell size
// EPSDEGEN = relative area below which a triangle is treated as degenerate

#define EPSPUSH 1.0e-6
#define EPSDEGEN 1.0e-12
#define BIGDSQ 1.0e300
#define DELTAPUSH 1024

#ifndef DEBUG_SPARTA_ENGULF
#define DEBUG_SPARTA_ENGULF 0
#endif
#ifndef DEBUG_SPARTA_WALLFRAME
#define DEBUG_SPARTA_WALLFRAME 0
#endif

// [ADDED] Phase 2: fix move/surf で移動中の表面グループ一覧
// used to look up the wall velocity of any moving surf element,
//   e.g. when 2 rotors are moved by 2 fix move/surf commands

static std::vector<MoveSurf *> moving_groups;

/* ---------------------------------------------------------------------- */

MoveSurf::MoveSurf(SPARTA *sparta) : Pointers(sparta)
{
  me = comm->me;
  nprocs = comm->nprocs;

  // pselect = 1 if point is moved, else 0

  if (domain->dimension == 2)
    memory->create(pselect,2*surf->nsurf,"move_surf:pselect");
  else
    memory->create(pselect,3*surf->nsurf,"move_surf:pselect");

  file = NULL;
  fp = NULL;

  // [ADDED] Phase 1: 壁内粒子の扱いと集計の初期化

  dim = domain->dimension;
  mode = 0;
  engulfflag = ENGULF_DELETE;
  wallvelflag = 0;
  for (int i = 0; i < NEGCOUNT; i++) eglast[i] = egtotal[i] = 0;
  warnflag = 0;
  npush_migrate = maxpush_migrate = 0;
  push_migrate = NULL;

  // [ADDED] Phase 2: 壁速度の初期化
  // fix move/surf creates this class on the step its motion starts

  frac_prev = 0.0;
  step_prev = update->ntimestep;
  velflag = registered = 0;
  omega = 0.0;
  vtrans[0] = vtrans[1] = vtrans[2] = 0.0;
  dmove = 0.0;
  rmaxaxis = -1.0;
}

/* ---------------------------------------------------------------------- */

MoveSurf::~MoveSurf()
{
  // [ADDED] Phase 2: 移動表面グループ一覧から外す

  if (registered) {
    for (int i = 0; i < (int) moving_groups.size(); i++)
      if (moving_groups[i] == this) {
        moving_groups.erase(moving_groups.begin() + i);
        break;
      }
  }

  memory->destroy(push_migrate);       // [ADDED] Phase 1
  memory->destroy(pselect);
  delete [] file;
  if (fp) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void MoveSurf::command(int narg, char **arg)
{
  if (!surf->exist)
    error->all(FLERR,"Cannot move_surf with no surf elements defined");

  if (surf->distributed)
    error->all(FLERR,
               "Cannot yet use move_surf with distributed surf elements");

  if (narg < 2) error->all(FLERR,"Illegal move_surf command");

  // process command-line args

  int igroup = surf->find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Move_surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  process_args(narg-1,&arg[1]);
  mode = 0;

  if (action == READFILE)
    error->all(FLERR,"Move_surf file option is not yet implemented");

  // perform surface move

  if (me == 0) {
    if (screen) fprintf(screen,"Moving surfs ...\n");
    if (logfile) fprintf(logfile,"Moving surfs ...\n");
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  dim = domain->dimension;

  // sort particles

  if (particle->exist) particle->sort();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // move line/tri points via chosen action by full amount

  if (dim == 2) move_lines(1.0,surf->lines);
  else move_tris(1.0,surf->tris);

  // remake list of surf elements I own
  // assign split cell particles to parent split cell
  // assign surfs to grid cells

  grid->unset_neighbors();
  grid->remove_ghosts();

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();
  grid->surf2grid(1);

  if (dim == 2) surf->check_point_near_surf_2d();
  else surf->check_point_near_surf_3d();

  if (dim == 2) surf->check_watertight_2d();
  else surf->check_watertight_3d();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // flag cells and corners as OUTSIDE or INSIDE
  // reallocate per grid cell arrays in per grid computes
  //   local grid cell counts could have changed due to split cell changes

  grid->set_inout();
  grid->type_check();

  // DEBUG
  //grid->debug();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  // remove particles as needed due to surface move

  bigint ndeleted;
  if (particle->exist) ndeleted = remove_particles();

  MPI_Barrier(world);
  double time6 = MPI_Wtime();

  double time_total = time6-time1;

  if (comm->me == 0) {
    if (screen) {
      if (particle->exist)
        fprintf(screen,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      // [ADDED] Phase 1: 退避した粒子数
      if (particle->exist && engulfflag == ENGULF_PROJECT)
        fprintf(screen,"  " BIGINT_FORMAT " particles pushed out of surfs\n",
                eglast[EG_PUSHED]);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  sort/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
    }
    if (logfile) {
      if (particle->exist)
        fprintf(logfile,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      // [ADDED] Phase 1: 退避した粒子数
      if (particle->exist && engulfflag == ENGULF_PROJECT)
        fprintf(logfile,"  " BIGINT_FORMAT " particles pushed out of surfs\n",
                eglast[EG_PUSHED]);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  sort/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   process command args for both move_surf and fix move/surf
------------------------------------------------------------------------- */

void MoveSurf::process_args(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal move surf command");

  int iarg = 0;
  if (strcmp(arg[0],"file") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal move surf command");
    action = READFILE;
    int n = strlen(arg[1]) + 1;
    file = new char[n];
    strcpy(file,arg[1]);
    n = strlen(arg[2]) + 1;
    entry = new char[n];
    strcpy(entry,arg[2]);
    iarg = 3;
  } else if (strcmp(arg[0],"trans") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal move surf command");
    action = TRANSLATE;
    delta[0] = input->numeric(FLERR,arg[1]);
    delta[1] = input->numeric(FLERR,arg[2]);
    delta[2] = input->numeric(FLERR,arg[3]);
    if (domain->dimension == 2 && delta[2] != 0.0)
      error->all(FLERR,"Invalid move surf translation for 2d simulation");
    iarg = 4;
  } else if (strcmp(arg[0],"rotate") == 0) {
    if (narg < 8) error->all(FLERR,"Illegal move surf command");
    action = ROTATE;
    theta = input->numeric(FLERR,arg[1]);
    rvec[0] = input->numeric(FLERR,arg[2]);
    rvec[1] = input->numeric(FLERR,arg[3]);
    rvec[2] = input->numeric(FLERR,arg[4]);
    origin[0] = input->numeric(FLERR,arg[5]);
    origin[1] = input->numeric(FLERR,arg[6]);
    origin[2] = input->numeric(FLERR,arg[7]);
    if (domain->dimension == 2 && (rvec[0] != 0.0 || rvec[1] != 0.0))
      error->all(FLERR,"Invalid move surf rotation for 2d simulation");
    if (rvec[0] == 0.0 && rvec[1] == 0.0 && rvec[2] == 0.0)
      error->all(FLERR,"Invalid move surf rotation");
    theta *= MY_PI/180.0;
    MathExtra::norm3(rvec);
    iarg = 8;
  } else error->all(FLERR,"Illegal move surf command");

  // optional args

  connectflag = 0;
  engulfflag = ENGULF_DELETE;          // [ADDED] Phase 1
  int wallvelset = -1;                 // [ADDED] Phase 2: -1 = not specified

  while (iarg < narg) {
    if (strcmp(arg[iarg],"connect") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal move surf command");
      if (strcmp(arg[iarg+1],"yes") == 0) connectflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) connectflag = 0;
      iarg += 2;

    // [ADDED] Phase 1: engulf delete|project
    // delete = legacy, delete all particles in cells with moved surfs
    // project = push particles inside surfs back into flow

    } else if (strcmp(arg[iarg],"engulf") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal move surf command");
      if (strcmp(arg[iarg+1],"delete") == 0) engulfflag = ENGULF_DELETE;
      else if (strcmp(arg[iarg+1],"project") == 0) engulfflag = ENGULF_PROJECT;
      else error->all(FLERR,"Illegal move surf command");
      iarg += 2;

    // [ADDED] Phase 2: wallvel yes|no
    // yes = all collisions with moved surfs use the wall velocity

    } else if (strcmp(arg[iarg],"wallvel") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal move surf command");
      if (strcmp(arg[iarg+1],"yes") == 0) wallvelset = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) wallvelset = 0;
      else error->all(FLERR,"Illegal move surf command");
      iarg += 2;
    } else error->all(FLERR,"Illegal move surf command");
  }

  // [ADDED] Phase 1: engulf project の適用条件

  if (engulfflag == ENGULF_PROJECT) {
    if (surf->implicit)
      error->all(FLERR,"Move surf engulf project requires explicit surfs");
    if (sparta->kokkos)
      error->all(FLERR,
                 "Move surf engulf project is not yet supported with Kokkos");
  }

  // [ADDED] Phase 2: 壁速度の適用条件と移動表面グループ一覧への登録
  // default = yes for engulf project, no for engulf delete
  // only fix move/surf (mode = 1) has a wall velocity
  // Update::setup() checks the list at the start of every run

  if (wallvelset < 0) wallvelflag = (engulfflag == ENGULF_PROJECT) ? 1 : 0;
  else wallvelflag = wallvelset;
  if (wallvelflag && engulfflag != ENGULF_PROJECT)
    error->all(FLERR,"Move surf wallvel yes requires engulf project");

  if (mode == 1 && wallvelflag && !registered) {
    moving_groups.push_back(this);
    registered = 1;
  }
}

/* ----------------------------------------------------------------------
   move points in lines via specified action
   each method sets pselect = 1 for moved points
   fraction = portion of full distance points should move
------------------------------------------------------------------------- */

void MoveSurf::move_lines(double fraction, Surf::Line *origlines)
{
  update_kinematics(fraction);         // [ADDED] Phase 2: 壁速度の更新

  if (connectflag && groupbit != 1) connect_2d_pre();

  if (action == READFILE) {
    readfile();
    update_points(fraction);
  }
  else if (action == TRANSLATE) translate_2d(fraction,origlines);
  else if (action == ROTATE) rotate_2d(fraction,origlines);

  if (connectflag && groupbit != 1) connect_2d_post();

  surf->compute_line_normal(0);

  // check that all points are still inside simulation box

  surf->check_point_inside(0);
}

/* ----------------------------------------------------------------------
   move points in triangles via specified action
   each method sets pselect = 1 for moved points
   fraction = portion of full distance points should move
------------------------------------------------------------------------- */

void MoveSurf::move_tris(double fraction, Surf::Tri *origtris)
{
  update_kinematics(fraction);         // [ADDED] Phase 2: 壁速度の更新

  if (connectflag && groupbit != 1) connect_3d_pre();

  if (action == READFILE) {
    readfile();
    update_points(fraction);
  }
  else if (action == TRANSLATE) translate_3d(fraction,origtris);
  else if (action == ROTATE) rotate_3d(fraction,origtris);

  if (connectflag && groupbit != 1) connect_3d_post();

  surf->compute_tri_normal(0);

  // check that all points are still inside simulation box

  surf->check_point_inside(0);
}

/* ----------------------------------------------------------------------
   read entry of new point coords from file
------------------------------------------------------------------------- */

void MoveSurf::readfile()
{
  /*
  int i;
  char line[MAXLINE];
  char *word,*eof;

  // open point file if necessary
  // if already open, will just continue scanning below
  // NOTE: allow for file name with wildcard char

  if (me == 0 && fp == NULL) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open move surf file %s",file);
      error->one(FLERR,str);
    }
  }

  // loop until section found that matches entry

  if (me == 0) {
    while (1) {
      if (fgets(line,MAXLINE,fp) == NULL)
        error->one(FLERR,"Did not find entry in move surf file");
      if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
      if (line[0] == '#') continue;                          // comment
      word = strtok(line," \t\n\r");
      if (strcmp(word,entry) != 0) continue;          // non-matching entry
      if (fgets(line,MAXLINE,fp) == NULL)
        error->one(FLERR,"Incompete entry in move surf file");
      word = strtok(line," \t\n\r");        // npoints value after entry
      nread = input->inumeric(FLERR,word);
    }
  }

  // allocate index and coord arrays for nread points

  MPI_Bcast(&nread,1,MPI_INT,0,world);

  if (oldcoord) {
    memory->destroy(readindex);
    memory->destroy(oldcoord);
    memory->destroy(newcoord);
    memory->create(readindex,nread,"move_surf:readindex");
    memory->create(oldcoord,nread,3,"move_surf:oldcoord");
    memory->create(newcoord,nread,3,"move_surf:newcoord");
  }

  // read nread point coords in entry
  // store old current point and new point coords so can move by fraction
  // rindex = ID (index) of this read-in point in master list
  // skip points that are out-of-range

  Surf::Point *pts = surf->pts;
  int npoint = surf->npoint;

  if (me == 0) {
    int id;
    double x,y,z;

    for (int i = 0; i < nread; i++) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one(FLERR,"Incomplete entry in move surf file");
      id = input->inumeric(FLERR,strtok(line," \t\n\r"));
      x = input->numeric(FLERR,strtok(NULL," \t\n\r"));
      y = input->numeric(FLERR,strtok(NULL," \t\n\r"));
      if (dim == 3) z = input->numeric(FLERR,strtok(NULL," \t\n\r"));
      else z = 0.0;
      if (id < 1 || id > npoint)
        error->one(FLERR,"Invalid point index in move surf file");
      id--;
      readindex[i] = id;
      oldcoord[i][0] = pts[id].x[0];
      oldcoord[i][1] = pts[id].x[1];
      oldcoord[i][2] = pts[id].x[2];
      newcoord[i][0] = x;
      newcoord[i][1] = y;
      newcoord[i][2] = z;
    }
  }

  // broadcast point info to all procs

  MPI_Bcast(readindex,nread,MPI_INT,0,world);
  MPI_Bcast(&oldcoord[0][0],3*nread,MPI_DOUBLE,0,world);
  MPI_Bcast(&newcoord[0][0],3*nread,MPI_DOUBLE,0,world);

  // pselect[I] = index of Ith surf point in nread points (for now)
  // NOTE: check that same surf point does not appear twice in nread list?

  for (i = 0; i < npoint; i++) pselect[i] = -1;
  for (i = 0; i < nread; i++) pselect[readindex[i]] = i;

  int *rflag;
  memory->create(rflag,nread,"move_surf:rflag");
  for (i =0; i < nread; i++) rflag[i] = 0;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nline = surf->nline;
  int ntri = surf->ntri;
  int p1,p2,p3;

  if (dim == 2) {
    for (i = 0; i < nline; i++) {
      if (!(lines[i].mask & groupbit)) continue;
      p1 = lines[i].p1;
      p2 = lines[i].p2;
      if (pselect[p1] >= 0) rflag[pselect[p1]] = 1;
      if (pselect[p2] >= 0) rflag[pselect[p2]] = 1;
    }
  } else {
    for (i = 0; i < ntri; i++) {
      if (!(tris[i].mask & groupbit)) continue;
      p1 = tris[i].p1;
      p2 = tris[i].p2;
      p3 = tris[i].p3;
      if (pselect[p1] >= 0) rflag[pselect[p1]] = 1;
      if (pselect[p2] >= 0) rflag[pselect[p2]] = 1;
      if (pselect[p3] >= 0) rflag[pselect[p3]] = 1;
    }
  }


  // pselect[I] = 1 if Ith surf point is moved by nread points, else 0

  for (i = 0; i < npoint; i++) pselect[i] = 0;
  for (i = 0; i < nread; i++)
    if (rflag[i]) pselect[readindex[i]] = 1;

  // clean up

  memory->destroy(rflag);
  */
}

/* ----------------------------------------------------------------------
   update points using info from file
------------------------------------------------------------------------- */

void MoveSurf::update_points(double fraction)
{
  /*
  int i;

  // update points by fraction of old to new move

  Surf::Point *pts = surf->pts;
  Surf::Point *p;

  for (i = 0; i < nread; i++) {
    if (pselect[readindex[i]] == 0) continue;
    p = &pts[readindex[i]];
    p->x[0] = oldcoord[i][0] + fraction * (newcoord[i][0]-oldcoord[i][0]);
    p->x[1] = oldcoord[i][1] + fraction * (newcoord[i][1]-oldcoord[i][1]);
    p->x[2] = oldcoord[i][2] + fraction * (newcoord[i][2]-oldcoord[i][2]);
  }
  */
}

/* ----------------------------------------------------------------------
   translate surf points in 2d
------------------------------------------------------------------------- */

void MoveSurf::translate_2d(double fraction, Surf::Line *origlines)
{
  double *p1,*p2,*op1,*op2;

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 2*nsurf; i++) pselect[i] = 0;

  double dx = fraction * delta[0];
  double dy = fraction * delta[1];

  for (int i = 0; i < nsurf; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    p1 = lines[i].p1;
    p2 = lines[i].p2;
    op1 = origlines[i].p1;
    op2 = origlines[i].p2;

    p1[0] = op1[0] + dx;
    p1[1] = op1[1] + dy;
    pselect[2*i] = 1;

    p2[0] = op2[0] + dx;
    p2[1] = op2[1] + dy;
    pselect[2*i+1] = 1;
  }
}

/* ----------------------------------------------------------------------
   translate surf points in 3d
------------------------------------------------------------------------- */

void MoveSurf::translate_3d(double fraction, Surf::Tri *origtris)
{
  double *p1,*p2,*p3,*op1,*op2,*op3;

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 3*nsurf; i++) pselect[i] = 0;

  double dx = fraction * delta[0];
  double dy = fraction * delta[1];
  double dz = fraction * delta[2];

  for (int i = 0; i < nsurf; i++) {
    if (!(tris[i].mask & groupbit)) continue;
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;
    op1 = origtris[i].p1;
    op2 = origtris[i].p2;
    op3 = origtris[i].p3;

    p1[0] = op1[0] + dx;
    p1[1] = op1[1] + dy;
    p1[2] = op1[2] + dz;
    pselect[3*i] = 1;

    p2[0] = op2[0] + dx;
    p2[1] = op2[1] + dy;
    p2[2] = op2[2] + dz;
    pselect[3*i+1] = 1;

    p3[0] = op3[0] + dx;
    p3[1] = op3[1] + dy;
    p3[2] = op3[2] + dz;
    pselect[3*i+2] = 1;
  }
}

/* ----------------------------------------------------------------------
   rotate surf points in 2d
------------------------------------------------------------------------- */

void MoveSurf::rotate_2d(double fraction, Surf::Line *origlines)
{
  double *p1,*p2,*op1,*op2;
  double q[4],d[3],dnew[3];
  double rotmat[3][3];

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 2*nsurf; i++) pselect[i] = 0;

  double angle = fraction * theta;
  MathExtra::axisangle_to_quat(rvec,angle,q);
  MathExtra::quat_to_mat(q,rotmat);

  for (int i = 0; i < nsurf; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    p1 = lines[i].p1;
    p2 = lines[i].p2;
    op1 = origlines[i].p1;
    op2 = origlines[i].p2;

    d[0] = op1[0] - origin[0];
    d[1] = op1[1] - origin[1];
    d[2] = op1[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p1[0] = dnew[0] + origin[0];
    p1[1] = dnew[1] + origin[1];
    pselect[2*i] = 1;

    d[0] = op2[0] - origin[0];
    d[1] = op2[1] - origin[1];
    d[2] = op2[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p2[0] = dnew[0] + origin[0];
    p2[1] = dnew[1] + origin[1];
    pselect[2*i+1] = 1;
  }
}

/* ----------------------------------------------------------------------
   rotate surf points in 3d
------------------------------------------------------------------------- */

void MoveSurf::rotate_3d(double fraction, Surf::Tri *origtris)
{
  double *p1,*p2,*p3,*op1,*op2,*op3;
  double q[4],d[3],dnew[3];
  double rotmat[3][3];

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 3*nsurf; i++) pselect[i] = 0;

  double angle = fraction * theta;
  MathExtra::axisangle_to_quat(rvec,angle,q);
  MathExtra::quat_to_mat(q,rotmat);

  for (int i = 0; i < nsurf; i++) {
    if (!(tris[i].mask & groupbit)) continue;
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;
    op1 = origtris[i].p1;
    op2 = origtris[i].p2;
    op3 = origtris[i].p3;

    d[0] = op1[0] - origin[0];
    d[1] = op1[1] - origin[1];
    d[2] = op1[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p1[0] = dnew[0] + origin[0];
    p1[1] = dnew[1] + origin[1];
    p1[2] = dnew[2] + origin[2];
    pselect[3*i] = 1;

    d[0] = op2[0] - origin[0];
    d[1] = op2[1] - origin[1];
    d[2] = op2[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p2[0] = dnew[0] + origin[0];
    p2[1] = dnew[1] + origin[1];
    p2[2] = dnew[2] + origin[2];
    pselect[3*i+1] = 1;

    d[0] = op3[0] - origin[0];
    d[1] = op3[1] - origin[1];
    d[2] = op3[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p3[0] = dnew[0] + origin[0];
    p3[1] = dnew[1] + origin[1];
    p3[2] = dnew[2] + origin[2];
    pselect[3*i+2] = 1;
  }
}

/* ----------------------------------------------------------------------
   add points in moved lines to hash
------------------------------------------------------------------------- */

void MoveSurf::connect_2d_pre()
{
  // hash for end points of moved lines
  // key = end point
  // value = global index (0 to 2*Nline-1) of the point
  // NOTE: could prealloc hash to correct size here

  hash = new MyHash();

  // add moved points to hash

  double *p1,*p2;
  OnePoint3d key;

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    p1 = lines[i].p1;
    p2 = lines[i].p2;
    key.pt[0] = p1[0]; key.pt[1] = p1[1]; key.pt[2] = 0.0;
    if (hash->find(key) == hash->end()) (*hash)[key] = 2*i+0;
    key.pt[0] = p2[0]; key.pt[1] = p2[1]; key.pt[2] = 0.0;
    if (hash->find(key) == hash->end()) (*hash)[key] = 2*i+1;
  }
}

/* ----------------------------------------------------------------------
   move points in lines connected to line points that were moved
------------------------------------------------------------------------- */

void MoveSurf::connect_2d_post()
{
  // check if non-moved points are in hash
  // if so, set their coords to matching point
  // set pselect for newly moved points so remove_particles() will work

  int m,value,j,jwhich;
  double *p[2],*q;
  OnePoint3d key;

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (lines[i].mask & groupbit) continue;
    p[0] = lines[i].p1;
    p[1] = lines[i].p2;

    for (m = 0; m < 2; m++) {
      key.pt[0] = p[m][0]; key.pt[1] = p[m][1]; key.pt[2] = 0.0;
      if (hash->find(key) != hash->end()) {
        value = (*hash)[key];
        j = value/2;
        jwhich = value % 2;
        if (jwhich == 0) q = lines[j].p1;
        else q = lines[j].p2;
        p[m][0] = q[0];
        p[m][1] = q[1];
        if (m == 0) pselect[2*i] = 1;
        else pselect[2*i+1] = 1;
      }
    }
  }

  // free the hash

  delete hash;
}

/* ----------------------------------------------------------------------
   add points in moved triangles to hash
------------------------------------------------------------------------- */

void MoveSurf::connect_3d_pre()
{
  // hash for corner points of moved triangles
  // key = corner point
  // value = global index (0 to 3*Ntri-1) of the point
  // NOTE: could prealloc hash to correct size here

  hash = new MyHash();

  // add moved points to hash

  double *p1,*p2,*p3;
  OnePoint3d key;

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (!(tris[i].mask & groupbit)) continue;
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;
    key.pt[0] = p1[0]; key.pt[1] = p1[1]; key.pt[2] = p1[2];
    if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+0;
    key.pt[0] = p2[0]; key.pt[1] = p2[1]; key.pt[2] = p2[2];
    if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+1;
    key.pt[0] = p3[0]; key.pt[1] = p3[1]; key.pt[2] = p3[2];
    if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+2;
  }
}

/* ----------------------------------------------------------------------
   move points in tris connected to tri points that were moved
------------------------------------------------------------------------- */

void MoveSurf::connect_3d_post()
{
  // check if non-moved points are in hash
  // if so, set their coords to matching point
  // set pselect for newly moved points so remove_particles() will work

  int m,value,j,jwhich;
  double *p[3],*q;
  OnePoint3d key;

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (tris[i].mask & groupbit) continue;
    p[0] = tris[i].p1;
    p[1] = tris[i].p2;
    p[2] = tris[i].p3;

    for (m = 0; m < 3; m++) {
      key.pt[0] = p[m][0]; key.pt[1] = p[m][1]; key.pt[2] = p[m][2];
      if (hash->find(key) != hash->end()) {
        value = (*hash)[key];
        j = value/3;
        jwhich = value % 3;
        if (jwhich == 0) q = tris[j].p1;
        else if (jwhich == 1) q = tris[j].p2;
        else q = tris[j].p3;
        p[m][0] = q[0];
        p[m][1] = q[1];
        p[m][2] = q[2];
        if (m == 0) pselect[3*i] = 1;
        else if (m == 1) pselect[3*i+1] = 1;
        else pselect[3*i+2] = 1;
      }
    }
  }

  // free the hash

  delete hash;
}

/* ----------------------------------------------------------------------
   remove particles in any cell that is now INSIDE or contains moved surfs
   surfs that moved determined by pselect for any of its points
   reassign particles in split cells to sub cell owner
   compress particles if any flagged for deletion
   NOTE: doc this logic better
------------------------------------------------------------------------- */

bigint MoveSurf::remove_particles()
{
  int isurf,nsurf;
  surfint *csurfs;

  dim = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;
  int delflag = 0;

  // [ADDED] Phase 1: 壁内粒子の集計（このプロセス分）

  bigint egcount[NEGCOUNT];
  for (int i = 0; i < NEGCOUNT; i++) egcount[i] = 0;
  npush_migrate = 0;

  for (int icell = 0; icell < nglocal; icell++) {

    // [ADDED] Phase 1: engulf project = 壁内粒子のみ幾何学的に退避
    // particles which cannot be rescued are only marked with icell = -1,
    //   they are removed by the single compress_rebalance() below
    // the legacy deletion below is used for engulf delete

    if (engulfflag == ENGULF_PROJECT) {
      if (push_out_cell(icell,egcount)) delflag = 1;
      continue;
    }

    // cell is inside surfs
    // remove particles in case it wasn't before

    if (cinfo[icell].type == INSIDE) {
      if (cinfo[icell].count) delflag = 1;
      particle->remove_all_from_cell(cinfo[icell].first);
      cinfo[icell].count = 0;
      cinfo[icell].first = -1;
      continue;
    }

    // cell has surfs or is split
    // if m < nsurf, loop over csurfs did not finish
    // which means cell contains a moved surf, so delete all its particles

    if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
      nsurf = cells[icell].nsurf;
      csurfs = cells[icell].csurfs;

      int m;
      if (dim == 2) {
        for (m = 0; m < nsurf; m++) {
          isurf = csurfs[m];
          if (pselect[2*isurf]) break;
          if (pselect[2*isurf+1]) break;
        }
      } else {
        for (m = 0; m < nsurf; m++) {
          isurf = csurfs[m];
          if (pselect[3*isurf]) break;
          if (pselect[3*isurf+1]) break;
          if (pselect[3*isurf+2]) break;
        }
      }

      if (m < nsurf) {
        if (cinfo[icell].count) delflag = 1;
        particle->remove_all_from_cell(cinfo[icell].first);
        cinfo[icell].count = 0;
        cinfo[icell].first = -1;
      }
    }

    if (cells[icell].nsplit > 1)
      grid->assign_split_cell_particles(icell);
  }

  // [ADDED] Phase 1: 他プロセスのセルへ退避した粒子の送受信
  // received particles which are not in the flow are marked for deletion

  if (engulfflag == ENGULF_PROJECT)
    if (receive_pushed(egcount)) delflag = 1;

  int nlocal_old = particle->nlocal;
  if (delflag) particle->compress_rebalance();
  bigint delta = nlocal_old - particle->nlocal;
  bigint ndeleted;
  MPI_Allreduce(&delta,&ndeleted,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // [ADDED] Phase 1: 退避で粒子の位置・所属セルが変わったので再ソートを要求
  // then sum tallies across procs

  if (engulfflag == ENGULF_PROJECT) {
    particle->sorted = 0;
    engulf_tally(egcount);
  }

  return ndeleted;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 幾何学的退避（1セル分）
   push particles in owned cell icell which are inside surfs back into flow
   only particles in INSIDE cells or cells with moved surfs are checked
   particles which cannot be rescued are marked with icell = -1,
     caller removes them with a single compress_rebalance()
   also assign remaining particles of a split cell to its sub cells
   return 1 if any particle was marked for deletion, else 0
------------------------------------------------------------------------- */

int MoveSurf::push_out_cell(int icell, bigint *egcount)
{
  int m,ip,isurf,inside;
  double xref[3];

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int *next = particle->next;
  int delflag = 0;

  // sub cells hold no particles, they are all in the parent split cell

  if (cells[icell].nsplit <= 0) return 0;

  // cell is entirely inside surfs, so every particle in it is inside
  // nearest surf is outside the cell, so particles always change cell

  if (cinfo[icell].type == INSIDE) {
    ip = cinfo[icell].first;
    while (ip >= 0) {
      if (!push_out_one(ip,icell,0,egcount)) delflag = 1;
      ip = next[ip];
    }
    cinfo[icell].count = 0;
    cinfo[icell].first = -1;
    return delflag;
  }

  // cell overlaps a moved non-transparent surf
  // test each particle geometrically, only particles inside surfs are moved
  // xref = point in flow, parity of surf crossings from particle to xref

  int nsurf = cells[icell].nsurf;
  if (nsurf && cinfo[icell].first >= 0) {
    surfint *csurfs = cells[icell].csurfs;
    int movedflag = 0;

    for (m = 0; m < nsurf; m++) {
      isurf = csurfs[m];
      if (dim == 2) {
        if (lines[isurf].transparent) continue;
        if (pselect[2*isurf] || pselect[2*isurf+1]) movedflag = 1;
      } else {
        if (tris[isurf].transparent) continue;
        if (pselect[3*isurf] || pselect[3*isurf+1] || pselect[3*isurf+2])
          movedflag = 1;
      }
      if (movedflag) break;
    }

    if (movedflag) {
      int refflag = grid->point_outside_surfs(icell,xref);
      ip = cinfo[icell].first;
      while (ip >= 0) {
        if (refflag)
          inside = particle_inside(icell,particle->particles[ip].x,xref);
        else inside = (cinfo[icell].volume <= 0.0);
        if (inside && !push_out_one(ip,icell,1,egcount)) delflag = 1;
        ip = next[ip];
      }
    }
  }

  // assign particles still in split cell to its sub cells
  // skip particles marked for deletion or moved to another cell

  if (cells[icell].nsplit > 1) {
    Particle::OnePart *particles = particle->particles;
    ip = cinfo[icell].first;
    while (ip >= 0) {
      if (particles[ip].icell == icell) {
        if (dim == 3) particles[ip].icell = update->split3d(icell,particles[ip].x);
        else particles[ip].icell = update->split2d(icell,particles[ip].x);
      }
      ip = next[ip];
    }
    cinfo[icell].count = 0;
    cinfo[icell].first = -1;
  }

  return delflag;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 1粒子の幾何学的退避
   particle ip in owned cell icell is inside a surf
   localflag = 1 if nearest surf may be among surfs of icell
   move it to its nearest surf point + small offset into the flow
   Phase 2: if approaching the moving wall, reflect it in the wall frame
   return 1 if particle survives
   return 0 if it cannot be rescued, it is then marked with icell = -1
------------------------------------------------------------------------- */

int MoveSurf::push_out_one(int ip, int icell, int localflag, bigint *egcount)
{
  int isurf,jcell;
  double dist;
  double xc[3],xnew[3],dir[3];

  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *p = &particle->particles[ip];

  egcount[EG_INSIDE]++;

  // nearest non-transparent surf element and closest point xc on it

  if (!nearest_surf(icell,p->x,localflag,xc,isurf,dist)) {
    p->icell = -1;
    egcount[EG_DEL_NOSURF]++;
    return 0;
  }

  // a particle engulfed by this move is at most dmove from the surf
  // allow 2 cells of slack, deeper particles are not moved

  double hcell = cell_size(icell);
  if (dist > dmove + 2.0*hcell) {
    p->icell = -1;
    egcount[EG_DEL_DEEP]++;
    return 0;
  }

  double *enorm;
  if (dim == 3) enorm = surf->tris[isurf].norm;
  else enorm = surf->lines[isurf].norm;

  // retreat direction:
  //   1st try = from particle thru closest point, valid at edges/corners
  //   2nd try = outward normal of the surf element
  // new position must be in a grid cell and outside all surfs
  // ghost cell = owned by another proc, particle is sent to it
  //   after all cells are processed and the owner checks it is in the flow

  double eps = EPSPUSH*hcell;
  int reason = EG_DEL_BLOCKED;
  jcell = -1;

  for (int itry = 0; itry < 2; itry++) {
    if (itry == 0) {
      if (dist <= eps) continue;
      dir[0] = (xc[0]-p->x[0]) / dist;
      dir[1] = (xc[1]-p->x[1]) / dist;
      dir[2] = (xc[2]-p->x[2]) / dist;
    } else {
      dir[0] = enorm[0];
      dir[1] = enorm[1];
      dir[2] = enorm[2];
    }
    if (dim == 2) dir[2] = 0.0;

    xnew[0] = xc[0] + eps*dir[0];
    xnew[1] = xc[1] + eps*dir[1];
    xnew[2] = xc[2] + eps*dir[2];
    if (dim == 2) xnew[2] = 0.0;

    jcell = locate_cell(icell,xnew);
    if (jcell < 0) {
      reason = EG_DEL_NOCELL;
      continue;
    }
    if (jcell >= grid->nlocal) break;
    if (point_in_flow(jcell,xnew)) break;
    reason = EG_DEL_BLOCKED;
    jcell = -1;
  }

  if (jcell < 0) {
    p->icell = -1;
    egcount[reason]++;
    return 0;
  }

  // Phase 2: collision with the moving wall at the closest point
  // only for fix move/surf, move_surf command has no wall velocity

  if (mode == 1) {
    if (!reflect_wall_frame(ip,isurf,xc,icell,egcount)) {
      particle->particles[ip].icell = -1;
      egcount[EG_DEL_DESTROYED]++;
      return 0;
    }
    p = &particle->particles[ip];
  }

  // commit new position and cell
  // ghost cell: keep ghost index, Comm::migrate_particles() sends particle
  // if particle stays in split cell icell, caller assigns its sub cell

  p->x[0] = xnew[0];
  p->x[1] = xnew[1];
  p->x[2] = xnew[2];

  if (jcell >= grid->nlocal) {
    if (npush_migrate == maxpush_migrate) {
      maxpush_migrate += DELTAPUSH;
      memory->grow(push_migrate,maxpush_migrate,"move_surf:push_migrate");
    }
    push_migrate[npush_migrate++] = ip;
    p->icell = jcell;
    egcount[EG_MIGRATED]++;
    return 1;
  }

  if (jcell != icell && cells[jcell].nsplit > 1) {
    if (dim == 3) jcell = update->split3d(jcell,xnew);
    else jcell = update->split2d(jcell,xnew);
  }
  p->icell = jcell;

  egcount[EG_PUSHED]++;
  return 1;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 最近傍の表面要素の探索
   find non-transparent surf element closest to point x in owned cell icell
   if localflag, search surfs of icell first
     result is certain if it is closer than any face of icell,
     since other surfs do not overlap icell
   else or if not certain, search all surfs (move_surf requires surfs
     not distributed, so every proc stores all of them)
   return 1 if found with xc = closest point, isurf = element, dist
   return 0 if no surf element exists
------------------------------------------------------------------------- */

int MoveSurf::nearest_surf(int icell, double *x, int localflag,
                           double *xc, int &isurf, double &dist)
{
  int m,i,j;
  double dsq,dbox;
  double xtmp[3];

  Grid::ChildCell *cells = grid->cells;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  double best = BIGDSQ;
  isurf = -1;

  if (localflag && cells[icell].nsurf) {
    int nsurf = cells[icell].nsurf;
    surfint *csurfs = cells[icell].csurfs;

    for (m = 0; m < nsurf; m++) {
      i = csurfs[m];
      if (dim == 2 && lines[i].transparent) continue;
      if (dim == 3 && tris[i].transparent) continue;
      dsq = closest_point_surf(i,x,xtmp);
      if (dsq < best) {
        best = dsq;
        isurf = i;
        xc[0] = xtmp[0]; xc[1] = xtmp[1]; xc[2] = xtmp[2];
      }
    }

    if (isurf >= 0) {
      double *lo = cells[icell].lo;
      double *hi = cells[icell].hi;
      dbox = BIGDSQ;
      for (j = 0; j < dim; j++) {
        if (x[j]-lo[j] < dbox) dbox = x[j]-lo[j];
        if (hi[j]-x[j] < dbox) dbox = hi[j]-x[j];
      }
      if (dbox > 0.0 && best <= dbox*dbox) {
        dist = sqrt(best);
        return 1;
      }
    }
  }

  // search all surf elements, skip those whose bounding box is too far

  int nsurf = surf->nsurf;

  for (i = 0; i < nsurf; i++) {
    if (dim == 2 && lines[i].transparent) continue;
    if (dim == 3 && tris[i].transparent) continue;
    if (bbox_distsq(i,x) >= best) continue;
    dsq = closest_point_surf(i,x,xtmp);
    if (dsq < best) {
      best = dsq;
      isurf = i;
      xc[0] = xtmp[0]; xc[1] = xtmp[1]; xc[2] = xtmp[2];
    }
  }

  if (isurf < 0) return 0;
  dist = sqrt(best);
  return 1;
}

/* ----------------------------------------------------------------------
   closest point C on line segment AB to point P
   return squared distance from P to C
------------------------------------------------------------------------- */

static double closest_point_segment(const double *a, const double *b,
                                    const double *p, double *c)
{
  double ab[3],ap[3],d[3];

  sub3(b,a,ab);
  sub3(p,a,ap);

  // degenerate segment of zero length: closest point is A

  double t = 0.0;
  double lensq = lensq3(ab);
  if (lensq > 0.0) {
    t = dot3(ap,ab) / lensq;
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;
  }

  c[0] = a[0] + t*ab[0];
  c[1] = a[1] + t*ab[1];
  c[2] = a[2] + t*ab[2];

  sub3(p,c,d);
  return lensq3(d);
}

/* ----------------------------------------------------------------------
   set C = A, return squared distance from P to C
------------------------------------------------------------------------- */

static double copy_point_distsq(const double *a, const double *p, double *c)
{
  double d[3];

  c[0] = a[0];
  c[1] = a[1];
  c[2] = a[2];

  sub3(p,c,d);
  return lensq3(d);
}

/* ----------------------------------------------------------------------
   closest point Q on triangle ABC to point P
   Voronoi region test, Ericson, Real-Time Collision Detection, sec 5.1.5
   every division is guarded, degenerate triangle uses its 3 edges
   return squared distance from P to Q
------------------------------------------------------------------------- */

static double closest_point_triangle(const double *a, const double *b,
                                     const double *c, const double *p,
                                     double *q)
{
  double ab[3],ac[3],ap[3],bp[3],cp[3],bc[3],d[3];
  double v,w,denom;

  sub3(b,a,ab);
  sub3(c,a,ac);
  sub3(p,a,ap);

  // vertex region of A

  double d1 = dot3(ab,ap);
  double d2 = dot3(ac,ap);
  if (d1 <= 0.0 && d2 <= 0.0) return copy_point_distsq(a,p,q);

  // vertex region of B

  sub3(p,b,bp);
  double d3 = dot3(ab,bp);
  double d4 = dot3(ac,bp);
  if (d3 >= 0.0 && d4 <= d3) return copy_point_distsq(b,p,q);

  // edge region of AB

  double vc = d1*d4 - d3*d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    denom = d1 - d3;
    v = (denom > 0.0) ? d1/denom : 0.0;
    q[0] = a[0] + v*ab[0];
    q[1] = a[1] + v*ab[1];
    q[2] = a[2] + v*ab[2];
    sub3(p,q,d);
    return lensq3(d);
  }

  // vertex region of C

  sub3(p,c,cp);
  double d5 = dot3(ab,cp);
  double d6 = dot3(ac,cp);
  if (d6 >= 0.0 && d5 <= d6) return copy_point_distsq(c,p,q);

  // edge region of AC

  double vb = d5*d2 - d1*d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    denom = d2 - d6;
    w = (denom > 0.0) ? d2/denom : 0.0;
    q[0] = a[0] + w*ac[0];
    q[1] = a[1] + w*ac[1];
    q[2] = a[2] + w*ac[2];
    sub3(p,q,d);
    return lensq3(d);
  }

  // edge region of BC

  double va = d3*d6 - d5*d4;
  if (va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0) {
    denom = (d4-d3) + (d5-d6);
    w = (denom > 0.0) ? (d4-d3)/denom : 0.0;
    sub3(c,b,bc);
    q[0] = b[0] + w*bc[0];
    q[1] = b[1] + w*bc[1];
    q[2] = b[2] + w*bc[2];
    sub3(p,q,d);
    return lensq3(d);
  }

  // face region
  // va+vb+vc = |AB x AC|^2, zero or tiny for a degenerate triangle

  denom = va + vb + vc;
  if (denom <= EPSDEGEN*lensq3(ab)*lensq3(ac)) {
    double qtmp[3];
    double best = closest_point_segment(a,b,p,q);
    double dsq = closest_point_segment(b,c,p,qtmp);
    if (dsq < best) {
      best = dsq;
      q[0] = qtmp[0]; q[1] = qtmp[1]; q[2] = qtmp[2];
    }
    dsq = closest_point_segment(c,a,p,qtmp);
    if (dsq < best) {
      best = dsq;
      q[0] = qtmp[0]; q[1] = qtmp[1]; q[2] = qtmp[2];
    }
    return best;
  }

  v = vb/denom;
  w = vc/denom;
  q[0] = a[0] + v*ab[0] + w*ac[0];
  q[1] = a[1] + v*ab[1] + w*ac[1];
  q[2] = a[2] + v*ab[2] + w*ac[2];
  sub3(p,q,d);
  return lensq3(d);
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 表面要素上の最近点
   closest point xc on surf element isurf to point x
   return squared distance
------------------------------------------------------------------------- */

double MoveSurf::closest_point_surf(int isurf, double *x, double *xc)
{
  if (dim == 3) {
    Surf::Tri *tri = &surf->tris[isurf];
    return closest_point_triangle(tri->p1,tri->p2,tri->p3,x,xc);
  }

  Surf::Line *line = &surf->lines[isurf];
  return closest_point_segment(line->p1,line->p2,x,xc);
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 表面要素の外接箱までの距離
   squared distance from point x to bounding box of surf element isurf
   lower bound on squared distance from x to the element
------------------------------------------------------------------------- */

double MoveSurf::bbox_distsq(int isurf, double *x)
{
  double lo,hi,delta;
  double *p1,*p2,*p3;

  if (dim == 3) {
    p1 = surf->tris[isurf].p1;
    p2 = surf->tris[isurf].p2;
    p3 = surf->tris[isurf].p3;
  } else {
    p1 = surf->lines[isurf].p1;
    p2 = surf->lines[isurf].p2;
    p3 = p2;
  }

  double dsq = 0.0;
  for (int j = 0; j < dim; j++) {
    lo = p1[j];
    if (p2[j] < lo) lo = p2[j];
    if (p3[j] < lo) lo = p3[j];
    hi = p1[j];
    if (p2[j] > hi) hi = p2[j];
    if (p3[j] > hi) hi = p3[j];
    if (x[j] < lo) delta = lo - x[j];
    else if (x[j] > hi) delta = x[j] - hi;
    else delta = 0.0;
    dsq += delta*delta;
  }

  return dsq;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 壁内判定（交差回数の偶奇）
   return 1 if point x in owned cell icell is inside surfs, else 0
   xref = point in flow within icell from Grid::point_outside_surfs()
   same test as Grid::outside_surfs(), but transparent surfs are skipped
     since particles pass thru them
------------------------------------------------------------------------- */

int MoveSurf::particle_inside(int icell, double *x, double *xref)
{
  int isurf,side;
  bool hitflag;
  double param;
  double xc[3];

  Grid::ChildCell *cells = grid->cells;
  surfint *csurfs = cells[icell].csurfs;
  int nsurf = cells[icell].nsurf;

  int cnt = 0;
  for (int m = 0; m < nsurf; m++) {
    isurf = csurfs[m];
    if (dim == 3) {
      Surf::Tri *tri = &surf->tris[isurf];
      if (tri->transparent) continue;
      hitflag = Geometry::
        line_tri_intersect_noeps(x,xref,tri->p1,tri->p2,tri->p3,
                                 tri->norm,xc,param,side);
    } else {
      Surf::Line *line = &surf->lines[isurf];
      if (line->transparent) continue;
      hitflag = Geometry::
        line_line_intersect(x,xref,line->p1,line->p2,line->norm,
                            xc,param,side);
    }
    if (hitflag) cnt++;
  }

  return cnt % 2;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 退避先が流体中か
   return 1 if point x in owned cell icell is in the flow, else 0
------------------------------------------------------------------------- */

int MoveSurf::point_in_flow(int icell, double *x)
{
  double xref[3];

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  if (cinfo[icell].type == INSIDE) return 0;
  if (cells[icell].nsurf == 0) return 1;

  // no flow point found means cell is only touched by surfs,
  //   it is then in the flow if it has flow volume

  if (!grid->point_outside_surfs(icell,xref))
    return (cinfo[icell].volume > 0.0);

  return !particle_inside(icell,x,xref);
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 退避先の所属セル
   return owned or ghost cell (unsplit or split, never a sub cell)
     containing x, ghost index is >= grid->nlocal
   icell = cell the particle is currently in, checked first
   return -1 if x is outside the box or beyond the ghost cells,
     caller deletes the particle
------------------------------------------------------------------------- */

int MoveSurf::locate_cell(int icell, double *x)
{
  int j;

  Grid::ChildCell *cells = grid->cells;
  double *lo = cells[icell].lo;
  double *hi = cells[icell].hi;

  int inflag = 1;
  for (j = 0; j < dim; j++)
    if (x[j] < lo[j] || x[j] > hi[j]) inflag = 0;
  if (inflag) return icell;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  for (j = 0; j < dim; j++)
    if (x[j] < boxlo[j] || x[j] > boxhi[j]) return -1;

  if (!grid->hashfilled) grid->rehash();

  return grid->id_find_child(0,0,boxlo,boxhi,x);
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 他プロセスのセルへ退避した粒子の送受信
   send particles pushed into ghost cells to the procs owning the cells
   Comm::migrate_particles() sets their icell to the owned cell index
   check each received particle is in the flow, else mark for deletion,
     then assign it to a sub cell if its cell is split
   return 1 if any particle was marked for deletion, else 0
------------------------------------------------------------------------- */

int MoveSurf::receive_pushed(bigint *egcount)
{
  // compress in Comm::migrate_particles() requires ascending indices

  std::sort(push_migrate,push_migrate+npush_migrate);
  int nstart = comm->migrate_particles(npush_migrate,push_migrate);
  npush_migrate = 0;

  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;
  int delflag = 0;

  for (int i = nstart; i < nlocal; i++) {
    int icell = particles[i].icell;
    if (!point_in_flow(icell,particles[i].x)) {
      particles[i].icell = -1;
      egcount[EG_DEL_BLOCKED]++;
      delflag = 1;
      continue;
    }
    if (cells[icell].nsplit > 1) {
      if (dim == 3) particles[i].icell = update->split3d(icell,particles[i].x);
      else particles[i].icell = update->split2d(icell,particles[i].x);
    }
    egcount[EG_PUSHED]++;
  }

  return delflag;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: セル寸法
   smallest side length of cell icell
------------------------------------------------------------------------- */

double MoveSurf::cell_size(int icell)
{
  double *lo = grid->cells[icell].lo;
  double *hi = grid->cells[icell].hi;

  double size = hi[0] - lo[0];
  if (hi[1]-lo[1] < size) size = hi[1] - lo[1];
  if (dim == 3 && hi[2]-lo[2] < size) size = hi[2] - lo[2];
  return size;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 1: 集計の全プロセス合計と警告
   sum tallies of this surf move across procs
   warn once if any particle could not be rescued
------------------------------------------------------------------------- */

void MoveSurf::engulf_tally(bigint *egcount)
{
  MPI_Allreduce(egcount,eglast,NEGCOUNT,MPI_SPARTA_BIGINT,MPI_SUM,world);
  for (int i = 0; i < NEGCOUNT; i++) egtotal[i] += eglast[i];

  bigint ndel = eglast[EG_DEL_NOSURF] + eglast[EG_DEL_DEEP] +
    eglast[EG_DEL_NOCELL] + eglast[EG_DEL_BLOCKED] + eglast[EG_DEL_DESTROYED];

  if (ndel && !warnflag) {
    warnflag = 1;
    if (me == 0) {
      char str[512];
      snprintf(str,512,"Move surf engulf project deleted " BIGINT_FORMAT
               " particles which could not be pushed out of surfs "
               "on step " BIGINT_FORMAT " (no surf " BIGINT_FORMAT
               ", too deep " BIGINT_FORMAT ", no cell " BIGINT_FORMAT
               ", blocked " BIGINT_FORMAT ", destroyed " BIGINT_FORMAT ")",
               ndel,update->ntimestep,eglast[EG_DEL_NOSURF],
               eglast[EG_DEL_DEEP],eglast[EG_DEL_NOCELL],
               eglast[EG_DEL_BLOCKED],eglast[EG_DEL_DESTROYED]);
      error->warning(FLERR,str);
    }
  }

#if DEBUG_SPARTA_ENGULF
  if (me == 0 && screen)
    fprintf(screen,"ENGULF step " BIGINT_FORMAT ": inside " BIGINT_FORMAT
            " pushed " BIGINT_FORMAT " via other procs " BIGINT_FORMAT
            " reflected " BIGINT_FORMAT " kept velocity " BIGINT_FORMAT
            " deleted " BIGINT_FORMAT "\n",
            update->ntimestep,eglast[EG_INSIDE],eglast[EG_PUSHED],
            eglast[EG_MIGRATED],eglast[EG_REFLECT],eglast[EG_KEEPVEL],ndel);
#endif
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 2: 壁速度の算出（真の速度 = 変位 / 経過時間）
   called at start of every surf move with fraction of full move
   fix move/surf: fraction = (step - step0) / Nlarge, so
     omega = theta * dfrac / (dstep * update->dt)
     vtrans = delta * dfrac / (dstep * update->dt)
   move_surf command moves instantaneously and has no wall velocity
   also set dmove = max displacement of any moved point in this move
------------------------------------------------------------------------- */

void MoveSurf::update_kinematics(double fraction)
{
  double dfrac;
  if (mode == 1) dfrac = fraction - frac_prev;
  else dfrac = fraction;

  // rigid rotation by angle phi moves a point at distance r from the axis
  //   by chord length 2 r sin(phi/2)

  if (action == TRANSLATE) dmove = fabs(dfrac) * len3(delta);
  else if (action == ROTATE) {
    if (rmaxaxis < 0.0) rmaxaxis = max_axis_distance();
    dmove = 2.0*rmaxaxis * fabs(sin(0.5*theta*dfrac));
  } else dmove = 0.0;

  if (mode != 1) return;

  // [REMOVED] Phase 2: 登録は process_args() へ移動（wallvel no のグループは登録しない）
  // if (!registered) {
  //   moving_groups.push_back(this);
  //   registered = 1;
  // }

  bigint dstep = update->ntimestep - step_prev;
  if (dstep > 0) {
    double dtime = dstep * update->dt;
    if (action == ROTATE) omega = theta*dfrac / dtime;
    else if (action == TRANSLATE) {
      vtrans[0] = delta[0]*dfrac / dtime;
      vtrans[1] = delta[1]*dfrac / dtime;
      vtrans[2] = delta[2]*dfrac / dtime;
    }
    if (action == ROTATE || action == TRANSLATE) velflag = 1;
  }

  frac_prev = fraction;
  step_prev = update->ntimestep;

#if DEBUG_SPARTA_WALLFRAME
  if (me == 0 && screen)
    fprintf(screen,"WALLFRAME step " BIGINT_FORMAT ": omega %g vtrans %g %g %g "
            "dmove %g\n",update->ntimestep,omega,
            vtrans[0],vtrans[1],vtrans[2],dmove);
#endif
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 2: 回転軸から最も遠い移動点の距離
------------------------------------------------------------------------- */

double MoveSurf::max_axis_distance()
{
  double d[3],dperp[3];
  double *pts[3];
  int npt;

  int dimension = domain->dimension;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  double rmax = 0.0;
  for (int i = 0; i < nsurf; i++) {
    if (dimension == 2) {
      if (!(lines[i].mask & groupbit)) continue;
      pts[0] = lines[i].p1;
      pts[1] = lines[i].p2;
      npt = 2;
    } else {
      if (!(tris[i].mask & groupbit)) continue;
      pts[0] = tris[i].p1;
      pts[1] = tris[i].p2;
      pts[2] = tris[i].p3;
      npt = 3;
    }
    for (int k = 0; k < npt; k++) {
      sub3(pts[k],origin,d);
      double dot = dot3(d,rvec);
      dperp[0] = d[0] - dot*rvec[0];
      dperp[1] = d[1] - dot*rvec[1];
      dperp[2] = d[2] - dot*rvec[2];
      double r = len3(dperp);
      if (r > rmax) rmax = r;
    }
  }

  return rmax;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 2: この表面グループの点 x での壁速度 u
------------------------------------------------------------------------- */

void MoveSurf::own_velocity(double *x, double *u)
{
  if (action == ROTATE) {
    double r[3],c[3];
    sub3(x,origin,r);
    cross3(rvec,r,c);
    u[0] = omega*c[0];
    u[1] = omega*c[1];
    u[2] = omega*c[2];
  } else if (action == TRANSLATE) {
    u[0] = vtrans[0];
    u[1] = vtrans[1];
    u[2] = vtrans[2];
  } else u[0] = u[1] = u[2] = 0.0;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 2: 表面要素 isurf 上の点 x での壁速度 u
   velocity of the fix move/surf whose group contains the element
   zero for static surfs or groups with wallvel no
------------------------------------------------------------------------- */

void MoveSurf::wall_velocity(int isurf, double *x, double *u)
{
  if (!wall_velocity_any(sparta,isurf,x,u)) u[0] = u[1] = u[2] = 0.0;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 2: 壁速度の参照（Update::move() 用の公開インタフェース）
   spa = SPARTA instance, isurf = surf element, x = point on it
   return 1 and set u = wall velocity at x if the element belongs to
     a fix move/surf group with wallvel yes and a known velocity
   return 0 if element is static, u is then not set
------------------------------------------------------------------------- */

int MoveSurf::wall_velocity_any(SPARTA *spa, int isurf, double *x, double *u)
{
  int mask;
  if (spa->domain->dimension == 3) mask = spa->surf->tris[isurf].mask;
  else mask = spa->surf->lines[isurf].mask;

  for (int i = 0; i < (int) moving_groups.size(); i++) {
    MoveSurf *ms = moving_groups[i];
    if (ms->sparta != spa || !ms->velflag) continue;
    if (!(mask & ms->groupbit)) continue;
    ms->own_velocity(x,u);
    return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 2: run 開始時の壁速度の有効判定（Update::setup() から）
   return 1 if any fix move/surf of this SPARTA instance uses wallvel yes
   error if one of its surfs uses surf_collide piston, which applies
     its own moving-wall reflection
   surf_collide translate/rotate keywords cannot be detected here,
     they must not be used on these surfs (wall velocity added twice)
------------------------------------------------------------------------- */

int MoveSurf::wallvel_setup(SPARTA *spa)
{
  int isc,mask;

  Surf *s = spa->surf;
  int dimension = spa->domain->dimension;
  int nsurf = s->nsurf;
  int flag = 0;

  for (int i = 0; i < (int) moving_groups.size(); i++) {
    MoveSurf *ms = moving_groups[i];
    if (ms->sparta != spa) continue;
    flag = 1;

    for (int j = 0; j < nsurf; j++) {
      if (dimension == 3) {
        mask = s->tris[j].mask;
        isc = s->tris[j].isc;
      } else {
        mask = s->lines[j].mask;
        isc = s->lines[j].isc;
      }
      if (!(mask & ms->groupbit) || isc < 0) continue;
      if (strcmp(s->sc[isc]->style,"piston") == 0)
        spa->error->all(FLERR,
                        "Move surf wallvel cannot be used with surf_collide piston");
    }
  }

  return flag;
}

/* ----------------------------------------------------------------------
   [ADDED] Phase 2: ガリレイ変換で挟んだ既存 collide() による反射
   particle ip hit surf element isurf at closest point xc
   w = v - u_w       transform to frame of moving wall
   w' = collide(w)   existing surface collision model, unchanged
   v' = w' + u_w     transform back to lab frame
   particle moving away from wall in wall frame keeps its velocity,
     it is inside the surf only due to the discrete surf move
   icell = owned cell the particle was in, used for surf tallies
   return 1 if particle survives, 0 if collision model destroyed it
------------------------------------------------------------------------- */

int MoveSurf::reflect_wall_frame(int ip, int isurf, double *xc, int icell,
                                 bigint *egcount)
{
  int isc;
  double *norm;
  double uw[3],w[3];

  if (dim == 3) {
    isc = surf->tris[isurf].isc;
    norm = surf->tris[isurf].norm;
  } else {
    isc = surf->lines[isurf].isc;
    norm = surf->lines[isurf].norm;
  }

  // no collision model, transparent model, or piston model:
  //   keep velocity, piston applies its own moving-wall reflection

  if (isc < 0) return 1;
  SurfCollide *sc = surf->sc[isc];
  if (sc->transparent || strcmp(sc->style,"piston") == 0) return 1;

  Particle::OnePart *p = &particle->particles[ip];

  wall_velocity(isurf,xc,uw);
  sub3(p->v,uw,w);

  if (dot3(w,norm) >= 0.0) {
    egcount[EG_KEEPVEL]++;
    return 1;
  }

  // lab-frame incident state for surface tallies

  Particle::OnePart iorig;
  int ntally = update->nsurf_tally;
  if (ntally) memcpy(&iorig,p,sizeof(Particle::OnePart));

  // w = v - u_w, collision at closest point
  // isr = -1, no surface chemistry for this delayed collision,
  //   so no collision model creates a particle

  p->x[0] = xc[0];
  p->x[1] = xc[1];
  p->x[2] = xc[2];
  p->v[0] = w[0];
  p->v[1] = w[1];
  p->v[2] = w[2];
  p->icell = icell;

  double dtremain = 0.0;
  int reaction = 0;
  Particle::OnePart *ipart = p;
  Particle::OnePart *jpart =
    sc->collide(ipart,dtremain,isurf,norm,-1,reaction);

  if (jpart)
    error->one(FLERR,"Move surf engulf project: "
               "surface collision model created a particle");

  // v' = w' + u_w

  if (ipart) {
    ipart->v[0] += uw[0];
    ipart->v[1] += uw[1];
    ipart->v[2] += uw[2];
  }

  if (ntally)
    for (int m = 0; m < ntally; m++)
      update->slist_active[m]->surf_tally(isurf,icell,reaction,
                                          &iorig,ipart,NULL);

  if (!ipart) return 0;

  egcount[EG_REFLECT]++;
  return 1;
}
