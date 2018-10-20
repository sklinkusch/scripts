/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/* Created by Anjuta version 1.2.3 */
/* Written by David Stevenson ( david at avoncliff dot com ) March 2006 */
/* Based on Daniele Mazzocchio's sudokiller.asm */
/* Compile with "gcc -o sudokiller sudokiller.c" */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define TRUE 1
#define FALSE 0

using namespace std;

void display_board(int* board);
int check(int cell, int num, int* board);
int guess(int cell, int* board);

int main(void){
 int* board = new int[81];
 cout << "Insert Sudoku, beginning with first row from left to right, space between numbers and '0' for unknown numbers\n";
 for(int i=0; i < 81; i++){
  cin >> board[i];
 }
	guess(0, board);
	display_board(board);
	return(0);
}

void display_board(int* board){
	int i, j;
	for(int k=0; k < 25; k++){
		if(k == 0 || k == 8 || k == 16 || k == 24){
			putchar(43);
		}else{
			putchar(45);
		}
	}
	putchar(10); putchar(13);
	for(i=0; i<9; i++){
		putchar(124); // prints vertical bar
		putchar(32); // prints space
		for(j=0; j<9; j++){ 
			putchar(0x30 + board[i*9+j]); // prints numbers from 1 to 9
			putchar(32); // prints space
			if(j%3 == 2){
				putchar(124); // prints vertical bar
				putchar(32); // prints space
			}
                }
		putchar(10); putchar(13);
		if(i%3 == 2){
			for(int k=0; k<25; k++){
				if(k == 0 || k == 8 || k == 16 || k == 24){
					putchar(43);
				}else{
					putchar(45);
				}
			}
			putchar(10); putchar(13);
		}
	}
}

int check(int cell, int num, int* board){
	int i, row, col;
	int block_row, block_col;
	row = 9 * (cell/9);
	col = cell%9;
	block_row = (row/27)*27;
	block_col = (col/3)*3;

	for(i=0; i<9; i++) if (board[row+i]==num) return FALSE;
	for(i=0; i<9; i++) if (board[col+(i*9)]==num) return FALSE;
	for(i=0; i<3; i++) if (board[block_row+block_col+i]==num) return FALSE;
	for(i=0; i<3; i++) if (board[(block_row+9)+block_col+i]==num) return FALSE;
	for(i=0; i<3; i++) if (board[(block_row+18)+block_col+i]==num) return FALSE;
	return TRUE;
}

int guess(int cell, int* board){
	int num;
	if (cell >= 81) return TRUE;
	if (board[cell]) return guess(cell+1, board);
	for (num=1; num<10; num++)
	{
		if (check(cell, num, board)){
			board[cell] = num;
			if (guess(cell+1, board)) return TRUE;
		}
	}
	board[cell] = 0;
	return FALSE;
}
