# Copyright (c) 2020-2023 Issaku Kanamori <kanamori-i@riken.jp>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.
#
# See the full license in the file "LICENSE".

CC = mpifccpx
#CPPFLAGS = 
#CFLAGS =
CFLAGS = 
#LDFLAGS = -ltofucom

SRC = rankmap_4d.c
OBJ = $(SRC:%.c=%.o)

OBJ1 = calc_rankid.o
OBJ2 = calc_rankid_reversed.o
PRG1 = rankmap_4d_lex
PRG2 = rankmap_4d_lex_reversed

PRG_GENERAL_1 = rankmap_4d_general_lex
PRG_GENERAL_2 = rankmap_4d_general_reversed
OBJ_GENERAL = rankmap_4d_general.o


all: $(PRG1) $(PRG2) $(PRG_GENERAL_1) $(PRG_GENERAL_2)

$(PRG1): $(OBJ) $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

$(PRG2): $(OBJ) $(OBJ2)
	$(CC) -o $@ $^ $(LDFLAGS)

$(PRG_GENERAL_1): $(OBJ_GENERAL) $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

$(PRG_GENERAL_2): $(OBJ_GENERAL) $(OBJ2)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o *.d *.lst

%.o: %.c
	$(CC) $(CFLAGS) -c -MMD $<
