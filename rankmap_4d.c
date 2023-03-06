/*
  4-dim rankmap generator for Fugaku
     Copyright(c) 2020, 2023, Issaku Kanamori <kanamori-i@riken.jp>

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 3
  of the License, or any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see http://www.gnu.org/licenses/.

  See the full license in the file "LICENSE".


    2020 Aug. 11 the first version
    2023 Mar.  6 added License description

 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <mpi-ext.h>
#include "config.h"

// global
int np;
int myrank;

typedef struct {
  int psize[4];
  int inner_dir;
} proc_dim;


void show_usage(char const * const *argv){
    printf("usage: %s P1 P2 P3 P4 [1234]\n", argv[0]);
    printf("       at least one of P1,P2,P3,P4 must be 4\n");
    printf("  ex. %s 8 4 4 4 4 --> 8x4x4x4 process lattice, 4th direction is the inner-node dirction\n", argv[0]);
    printf("  ex. %s 8 4 4 4   --> 8x4x4x4 process lattice, 2nd (1st \"4\") is the inner-node dirction\n", argv[0]);
}

// defined in calc_rankid.c
int calc_rankid(const int *coords, const int *psize);
extern const char* rankmap_name;

/**************************************************

  utility functions

**************************************************/
void safe_abort_(int status, const char* file, int line){
  if(myrank==0){
    printf("safe_abort is called in %s, at line %d:  status=%d\n", file, line, status);
  }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(status);
}

#define safe_abort(status) safe_abort_(status, __FILE__, __LINE__);

void check_error(const int rc, const int success, const char* msg){
  int flag=0;
  int recv=0;
  if(rc != success){
    flag=1;
  }
  MPI_Allreduce(&flag, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(recv>0){
    fprintf(stderr, "error at %s\n", msg);
    safe_abort(EXIT_FAILURE);
  }
}


void get_param(proc_dim *dim, const int argc, char const * const *argv){

  int *proc=dim->psize;  // alias
  proc[0]=atoi(argv[1]);  // px
  proc[1]=atoi(argv[2]);  // py
  proc[2]=atoi(argv[3]);  // pz
  proc[3]=atoi(argv[4]);  // pt

  if(np != proc[0]*proc[1]*proc[2]*proc[3] ){
    if(myrank==0){
      printf("np=%d != p1 x p2 x p3 x p4\n", np);
      printf("p1,p2,p3,p4=%d,%d,%d,%d\n",proc[0], proc[1], proc[2], proc[3]);
    }
    safe_abort(EXIT_FAILURE);
  }

  int dir=-1;
  if(argc>5){
    char c=argv[5][0];
    switch(c){
    case '1':
      dir=0;
      break;
    case '2':
      dir=1;
      break;
    case '3':
      dir=2;
      break;
    case '4':
      dir=3;
      break;
    default:
      if(myrank==0){
        printf("bad dir char: must be [1234] but given as %c\n", c);
      }
      safe_abort(EXIT_FAILURE);
    }
    if(proc[dir] !=4){
      if(myrank==0){
        printf("dir char is %c but proc[%d] = %d != 4\n", c, dir, proc[dir]);
      }
      safe_abort(EXIT_FAILURE);
    }
  } else {
    for(int i=0; i<4; i++){
      if(proc[i]==4) {
        dir=i;
        break;
      }
    }
    if(dir<0){
      if(myrank==0){
        printf("none of the proc size is 4 (at least one must be 4): %d,%d,%d,%d\n",
               proc[0],proc[1],proc[2],proc[3]);
      }
      safe_abort(EXIT_FAILURE);
    }
  }
  assert(dir == 0 || dir == 1 || dir == 2 || dir == 3);
  dim->inner_dir=dir;
  return;
}


void set_direction_map(int *dirmap, const proc_dim *dim, const int *shape_fjmpi){
  // dirmap[dir]      (dir: 0--3) 
  //   = -1   if dir is the in-node direction
  //  or
  //   = ( direction in the given 3-dim topology)
  int flag[4]={0};
  for(int i=0; i<4; i++){
    if(i==dim->inner_dir){
      dirmap[i]=-1;
      flag[3]++;
      continue;
    }
    for(int fjdir=0; fjdir<3; fjdir++){
      if(flag[fjdir]>0){ continue; }
      if(dim->psize[i] == shape_fjmpi[fjdir]){
        dirmap[i]=fjdir;
        flag[fjdir]++;
        break;
      }
    }
  }

  // sanity check
  if(flag[0]*flag[1]*flag[2]*flag[3] != 1){

    if(myrank==0){
      fprintf(stderr, "something is wrong in the process size, cannot map the process to the given topology\n");
      fprintf(stderr, " required 4-dim process size: %d %d %d %d\n",
              dim->psize[0], dim->psize[1], dim->psize[2], dim->psize[3]);
      fprintf(stderr, " in-node direction [1-4]: %d\n", dim->inner_dir+1);
      fprintf(stderr, " 3-dim node shape: %d %d %d\n", shape_fjmpi[0], shape_fjmpi[1], shape_fjmpi[2]);
    }
    safe_abort(EXIT_FAILURE);
  }
}


/********************************************************
 * actual work
 *   obtain 3 dim MPI coodinate and map to 1 dim rank id
 *   for a suitable 4 dim map
 *   N.B. the map from 4dim to 1dim is defeind in calc_rankid()
 *
 ********************************************************/
void set_rankmap(int *rank_list, const proc_dim *dim){
  int list_size=3*np;

  // clear
  for(int i=0; i<list_size; i++){
    rank_list[i]=-1;
  }

  // obatin the 3dim rank coordinate
  int rc;
  int dim_fjmpi;
  rc = FJMPI_Topology_get_dimension(&dim_fjmpi);
  check_error(rc, FJMPI_SUCCESS, "FJMPI_Toplogy_get_dimension");
  if(dim_fjmpi != 3){
    if(myrank == 0){
      printf("rank topolog must be 3-dim but given dimension is %d\n", dim_fjmpi);
    }
    safe_abort(EXIT_FAILURE);
  }

  int coords_fjmpi[3];
  rc = FJMPI_Topology_get_coords(MPI_COMM_WORLD, myrank, FJMPI_LOGICAL, dim_fjmpi, coords_fjmpi);
  check_error(rc, FJMPI_SUCCESS, "FJMPI_Toplogy_get_coords");
  int shape_fjmpi[3];
  rc = FJMPI_Topology_get_shape(shape_fjmpi, shape_fjmpi+1, shape_fjmpi+2);
  check_error(rc, FJMPI_SUCCESS, "FJMPI_Toplogy_get_shape");
  if(myrank==0){
    printf("shape of FJMPI: %d %d %d\n", shape_fjmpi[0], shape_fjmpi[1], shape_fjmpi[2]);
    printf("using rankmap: %s\n", rankmap_name);
  }

  int dirmap[4];
  set_direction_map(dirmap, dim, shape_fjmpi);

  //  int dim3=0;
  int coords[4];
  for(int i=0; i<4; i++){
    if(dirmap[i]<0){
      coords[i] = myrank % 4;
    } else {
      coords[i] = coords_fjmpi[dirmap[i]];
    }
  }

  int rankid=calc_rankid(coords, dim->psize);
  //  printf("I am %d: rankid=%d, coords=%d,%d,%d,%d\n", myrank, rankid, coords[0], coords[1], coords[2], coords[3]);
  //  fflush(0);

  // prepare the work area and set the coordiante of this process
  int *work=malloc(sizeof(int)*list_size);
  for(int i=0; i<list_size; i++){
    work[i]=0;
  }
  int offset=3*rankid;
  for(int i=0; i<dim_fjmpi; i++){
    work[offset+i]=coords_fjmpi[i];
  }

  // obtain the coordinate of all the process
  MPI_Allreduce(work, rank_list, list_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  free(work);

  // sanity check
  for(int i=0; i<np; i++){
    for(int dim3=0; dim3<3; dim3++){
      int r=rank_list[3*i+dim3];
      if(r < 0 || r >=shape_fjmpi[dim3]){
        if(myrank==0){
          fprintf(stderr, "cannot happen: i=%d, rank_list[3*i + %d]=%d\n", i, dim3,r);
        }
        safe_abort(EXIT_FAILURE);
      }
    } // dim3
  } // i

  return;
}


/***********************************************************
 * out put the rankmap to a file
 *   the output filename is defined with macro
 ***********************************************************/
void output_rankmap(const int *rank_list, const proc_dim *dim){

  FILE *fp;
  int err=0;
  const char *filename=RANK_MAP_FILE;
  if(myrank==0){
    printf("rank map file: %s\n", filename);
    fp=fopen(filename, "w");
    if(!fp){
      err=1;
      fprintf(stderr, "cannot open the output file: %s\n", filename);
    }
  }
  check_error(err, 0, "openning the output file");
  if(myrank==0){
    for(int i=0; i<np; i++){
      fprintf(fp,"(%d,%d,%d)\n", rank_list[3*i], rank_list[3*i+1], rank_list[3*i+2]);
    }
    fclose(fp);
  }
  return;
}


int main(int argc, char** argv){
  // initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  // read parameters
  if(argc<5){
    if(myrank==0){
      show_usage(argv);
    }
    safe_abort(EXIT_FAILURE);
  }
  proc_dim proc;
  get_param(&proc, argc, argv);

  // allocate rankmap list
  int *rank_list=malloc(sizeof(int)*3*np);

  // generate rankmap
  set_rankmap(rank_list, &proc);

  // output the rankmap to file
  output_rankmap(rank_list, &proc);

  // reallocate
  free(rank_list);

  // done
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("finished: rankmap_4d.\n");
    fflush(stdout);
  }
  MPI_Finalize();
  return 0;
}
