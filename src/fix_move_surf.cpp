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

#include "string.h"
#include <cstdio>
#include "fix_move_surf.h"
#include "move_surf.h"
#include "update.h"
#include "comm.h"
#include "input.h"
#include "grid.h"
#include "surf.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

/* ---------------------------------------------------------------------- */

FixMoveSurf::FixMoveSurf(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix move/surf command");

  if (!surf->exist)
    error->all(FLERR,"Cannot fix move/surf with no surf elements are defined");

  if (surf->distributed)
    error->all(FLERR,
               "Cannot yet use fix move/surf with distributed surf elements");

  scalar_flag = 1;
  global_freq = 1;

  // create instance of MoveSurf class

  me = comm->me;
  nprocs = comm->nprocs;

  movesurf = new MoveSurf(sparta);
  movesurf->mode = 1;

  // parse and check arguments using MoveSurf class

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute surf group ID does not exist");
  movesurf->groupbit = surf->bitmask[igroup];

  nevery = input->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix move/surf command");
  if (update->ntimestep % nevery)
    error->all(FLERR,"Current timestep must be multiple of "
               "fix move/surf nevery");

  nlarge = input->inumeric(FLERR,arg[4]);
  if (nlarge < 0) error->all(FLERR,"Illegal fix move/surf command");
  if (nlarge % nevery)
    error->all(FLERR,"Fix move/surf nlarge must be multiple of nevery");

  movesurf->process_args(narg-5,&arg[5]);

  dim = domain->dimension;
  ntimestep_original = update->ntimestep;

  // make copy of original lines or tris
  // to pass to movesurf->move_lines() or movesurf->move_tris()

  nsurf = surf->nsurf;
  origlines = NULL;
  origtris = NULL;

  if (dim == 2) {
    origlines = (Surf::Line *)
      memory->smalloc(nsurf*sizeof(Surf::Line),"fix/move/surf:origlines");
    memcpy(origlines,surf->lines,nsurf*sizeof(Surf::Line));
  } else if (dim == 3) {
    origtris = (Surf::Tri *)
      memory->smalloc(nsurf*sizeof(Surf::Tri),"fix/move/surf:origtris");
    memcpy(origtris,surf->tris,nsurf*sizeof(Surf::Tri));
  }

  // initial output

  ndeleted = 0;
}

/* ---------------------------------------------------------------------- */

FixMoveSurf::~FixMoveSurf()
{
  delete movesurf;
  memory->sfree(origlines);
  memory->sfree(origtris);
}

/* ---------------------------------------------------------------------- */

int FixMoveSurf::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMoveSurf::init()
{
  if (nsurf != surf->nsurf)
    error->all(FLERR,"Number of surface elements changed in fix move/surf");

  // NOTE: first read of file ?
  //       what about on successive run
  //       for both file and trans/rotate
}

/* ----------------------------------------------------------------------
   check for zero-volume cells before collide()
   this prevents "Collision cell volume is zero" errors
------------------------------------------------------------------------- */

void FixMoveSurf::start_of_step()
{
  // Remove particles in zero-volume cells before collide()
  // This is critical to prevent collision errors
  // Note: This is called before particle->sort() in update.cpp,
  // but we also need to check after sort() since sort() may assign
  // particles to zero-volume cells. However, we can't check after
  // sort() here, so we rely on the fact that sort() should not
  // assign particles to zero-volume cells if the grid is properly set up.
  // If errors still occur, we may need to modify update.cpp to check
  // after sort() but before collide().
  if (particle->exist) {
    Particle::OnePart *particles = particle->particles;
    Grid::ChildInfo *cinfo = grid->cinfo;
    int nlocal = particle->nlocal;
    int delflag = 0;
    
    for (int i = 0; i < nlocal; i++) {
      int icell = particles[i].icell;
      // CRITICAL: Check against maxlocal (array size) not nlocal
      if (icell >= 0 && icell < grid->maxlocal) {
        if (icell < grid->nlocal && cinfo[icell].volume == 0.0) {
          particles[i].icell = -1;
          delflag = 1;
        } else if (icell >= grid->nlocal) {
          // icell is within array bounds but exceeds nlocal
          // This can happen if grid was rebalanced - mark for deletion
          particles[i].icell = -1;
          delflag = 1;
        }
      } else if (icell >= grid->maxlocal) {
        // icell exceeds array bounds - mark for deletion
        particles[i].icell = -1;
        delflag = 1;
      }
    }
    
    if (delflag) {
      particle->compress_rebalance();
    }
  }
}

/* ----------------------------------------------------------------------
   move surface points incrementally
------------------------------------------------------------------------- */

void FixMoveSurf::end_of_step()
{
  // fraction = fraction of Nlarge represented by this timestep

  bigint nelapsed = update->ntimestep - ntimestep_original;
  double fraction = 1.0*nelapsed / nlarge;

  // movesurf moves surface vertices

  if (dim == 2) movesurf->move_lines(fraction,origlines);
  else movesurf->move_tris(fraction,origtris);

  // assign split cell particles to parent split cell
  // assign surfs to grid cells

  grid->unset_neighbors();
  grid->remove_ghosts();

  if (grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();

  grid->surf2grid(1,0);

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  // CRITICAL: Rebuild hash table after grid changes
  // This ensures hash table contains valid cell indices after grid rebalancing
  if (grid->hashfilled) {
    grid->rehash();
  }
  
  // CRITICAL: After grid rebalancing, validate and fix all particle icell values
  // Grid rebalancing may change maxlocal, making existing icell values invalid
  if (particle->exist) {
    Particle::OnePart *particles = particle->particles;
    int nlocal = particle->nlocal;
    int maxlocal = grid->maxlocal;
    int delflag = 0;
    
    for (int i = 0; i < nlocal; i++) {
      int icell = particles[i].icell;
      // Mark particles with invalid icell for deletion
      // icell must be < maxlocal (array bounds) and < nlocal (owned cells only)
      if (icell < 0 || icell >= maxlocal || icell >= grid->nlocal) {
        particles[i].icell = -1;
        delflag = 1;
      }
    }
    
    if (delflag) {
      particle->compress_rebalance();
      // Re-fetch particles pointer after compress_rebalance()
      particles = particle->particles;
      nlocal = particle->nlocal;
    }
  }
  
  grid->reset_neighbors();
  comm->reset_neighbors();

  // flag cells and corners as OUTSIDE or INSIDE

  grid->set_inout();
  grid->type_check(0);

  // 1. Detect trapped particles (Logging only)
  int *trapped_list = NULL;
  int ntrapped = 0;
  if (particle->exist) {
    movesurf->detect_trapped_particles(trapped_list, ntrapped);
  }

  // 2. Remove particles (NO PUSH OUT)
  // ここで push_out_particles を呼んでいた行を完全に削除すること！
  if (particle->exist) {
    ndeleted = movesurf->remove_particles();
    if (comm->me == 0) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
    }
    // CRITICAL: sort particles after remove_particles()
    // remove_particles() calls assign_split_cell_particles() which changes icell
    // but does not update the linked list, so we need to sort
    particle->sort();
  } else {
    ndeleted = 0;
  }

  // 3. Free trapped_list if allocated
  if (trapped_list) memory->sfree(trapped_list);

  // notify all classes that store per-grid data that grid may have changed

  grid->notify_changed();
}

/* ----------------------------------------------------------------------
   return particles deleted on last surface move
------------------------------------------------------------------------- */

double FixMoveSurf::compute_scalar()
{
  return (double) ndeleted;
}
