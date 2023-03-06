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

// change this function for different rankmap
int calc_rankid(const int *coords, const int *psize){
  return coords[0] + psize[0]*(coords[1] + psize[1]*(coords[2] + psize[2]*coords[3]));
}

void get_rank_coord(int *coords, int rank, const int *psize){
  int tmp=rank;
  for(int i=0; i<3; i++){
    coords[i]=tmp % psize[i];
    tmp /= psize[i];
  }
  coords[3]=tmp;
}

// for output log
const char* rankmap_name="lexical rankmap";
