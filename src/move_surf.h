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

#ifdef COMMAND_CLASS

CommandStyle(move_surf,MoveSurf)

#else

#ifndef SPARTA_MOVE_SURF_H
#define SPARTA_MOVE_SURF_H

#include "stdio.h"
#include "pointers.h"
#include "surf.h"
#include "hash3.h"

namespace SPARTA_NS {

class MoveSurf : protected Pointers {
 public:
  int mode;                 // 0 = move_surf command, 1 = fix move/surf command
  int groupbit;             // FixMoveSurf sets surf group

  // [ADDED] Phase 1: 壁内粒子の扱い（engulf キーワード）
  int engulfflag;           // 0 = delete all particles in cells with moved surfs
                            // 1 = push particles inside surfs back into flow

  // hash for surface element points

#include "hash_options.h"

#ifdef SPARTA_MAP
  typedef std::map<OnePoint3d,int> MyHash;
#elif SPARTA_UNORDERED_MAP
  typedef std::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
#else
  typedef std::tr1::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
#endif

  MyHash *hash;

  MoveSurf(class SPARTA *);
  ~MoveSurf();
  void command(int, char **);
  void process_args(int, char **);
  void move_lines(double, Surf::Line *);
  void move_tris(double, Surf::Tri *);
  bigint remove_particles();

  // [ADDED] Phase 2: 移動壁の壁速度（Update::move() の衝突から参照）

  static int wallvel_setup(class SPARTA *);
  static int wall_velocity_any(class SPARTA *, int, double *, double *);

 private:
  int me,nprocs;
  int dim,action,connectflag;
  char *file,*entry;
  double theta;
  double delta[3],rvec[3],origin[3];
  FILE *fp;

  int *pselect;                    // 1 if point is moved, else 0

  int nread;
  int *readindex;
  double **oldcoord,**newcoord;

  // [ADDED] Phase 1: 壁内粒子の集計（全プロセス合計値）

  enum{EG_INSIDE,EG_PUSHED,EG_MIGRATED,EG_REFLECT,EG_KEEPVEL,
       EG_DEL_NOSURF,EG_DEL_DEEP,EG_DEL_NOCELL,EG_DEL_BLOCKED,
       EG_DEL_DESTROYED,NEGCOUNT};

  bigint eglast[NEGCOUNT];         // tallies for last surf move
  bigint egtotal[NEGCOUNT];        // cumulative tallies
  int warnflag;                    // 1 if deletion warning already issued

  int npush_migrate;               // # of particles pushed into cells of
  int maxpush_migrate;             //   other procs, sent after all cells
  int *push_migrate;               //   are processed

  // [ADDED] Phase 2: 移動表面の剛体運動学（壁速度）

  int wallvelflag;                 // 1 if all collisions with this group
                                   //   are done in the frame of the moving
                                   //   wall, its surf_collide models must
                                   //   not use translate/rotate keywords
  double frac_prev;                // fraction of full move at previous move
  bigint step_prev;                // timestep of previous move
  int velflag;                     // 1 once wall velocity is known
  int registered;                  // 1 if in list of moving surf groups
  double omega;                    // angular velocity around rvec (rad/time)
  double vtrans[3];                // translational velocity (distance/time)
  double dmove;                    // max displacement of any point in last move
  double rmaxaxis;                 // max distance of group points from axis

  void readfile();
  void update_points(double);
  void translate_2d(double, Surf::Line *);
  void translate_3d(double, Surf::Tri *);
  void rotate_2d(double, Surf::Line *);
  void rotate_3d(double, Surf::Tri *);
  void connect_2d_pre();
  void connect_2d_post();
  void connect_3d_pre();
  void connect_3d_post();

  // [ADDED] Phase 1: 幾何学的退避と遅延削除

  int push_out_cell(int, bigint *);
  int push_out_one(int, int, int, bigint *);
  int nearest_surf(int, double *, int, double *, int &, double &);
  double closest_point_surf(int, double *, double *);
  double bbox_distsq(int, double *);
  int particle_inside(int, double *, double *);
  int point_in_flow(int, double *);
  int locate_cell(int, double *);
  double cell_size(int);
  int receive_pushed(bigint *);
  void engulf_tally(bigint *);

  // [ADDED] Phase 2: ガリレイ変換による壁速度の伝達

  void update_kinematics(double);
  double max_axis_distance();
  void own_velocity(double *, double *);
  void wall_velocity(int, double *, double *);
  int reflect_wall_frame(int, int, double *, int, bigint *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Move surf engulf project requires explicit surfs

The push-out of particles inside surfs uses explicit line/tri elements.

E: Move surf engulf project is not yet supported with Kokkos

The push-out calls surface collision models on the host.

W: Move surf engulf project deleted particles which could not be pushed out of surfs

Some particles inside surfs had no surf element, were deeper than the
surf move distance, or could not be placed in the flow of any grid cell.
They were deleted and are counted in the fix move/surf output.

E: Move surf wallvel yes requires engulf project

Collisions in the frame of the moving wall are only consistent if
particles inside surfs are pushed out instead of deleted.

E: Move surf wallvel cannot be used with surf_collide piston

The piston model applies its own moving-wall reflection.

E: Move surf engulf project: surface collision model created a particle

The delayed collision of a pushed-out particle is done without surface
chemistry, so no collision model should create a particle.

*/
