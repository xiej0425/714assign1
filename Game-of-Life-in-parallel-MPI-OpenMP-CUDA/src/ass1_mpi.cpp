#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <vector>
#include <iostream>

using namespace std;

int main (int argc, char* argv[]){
    int mpirank, mpisize;

    // Initialize the MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    int mpiroot = 0;
    int nRows, nCols, nTime;

    if (mpirank == mpiroot) {
        if (argc != 4){
            cerr << "incorrect number of input, expected 3" << endl;
            exit(1);
        }
        nRows = atoi(argv[1]);
        nCols = atoi(argv[2]);
        nTime = atoi(argv[3]);
        if (nRows < 1 || nCols < 1 || nTime < 1) {
            cerr << "incorrect input range" << endl;
            exit(1);
        }
    }

    MPI_Bcast(&nRows, 1, MPI_INT, mpiroot, MPI_COMM_WORLD);
    MPI_Bcast(&nCols, 1, MPI_INT, mpiroot, MPI_COMM_WORLD);
    MPI_Bcast(&nTime, 1, MPI_INT, mpiroot, MPI_COMM_WORLD);

    int nRowsLocal = nRows / mpisize;
    if (mpirank == mpisize - 1){
        nRowsLocal += nRows % mpisize;
    }

    int nRowsLocalWithGhost = nRowsLocal + 2;
    int nColsWithGhost = nCols + 2;

    vector<vector<int> > currGrid(nRowsLocalWithGhost, vector<int>
            (nColsWithGhost, 0));
    vector<vector<int> > nextGrid(nRowsLocalWithGhost, vector<int>
            (nColsWithGhost, 0));

    for (int iRow = 1; iRow <= nRows; ++iRow){
        for (int iCol = 1; iCol <= nCols; ++iCol){
            currGrid[iRow][iCol] = rand() % 2;
        }
    }

    int upperNeighbor = (mpirank == 0) ? mpisize - 1 : mpirank - 1;
    int lowerNeighbor = (mpirank == mpisize - 1) ? 0 : mpirank + 1;

    const int ALIVE = 1;
    const int DEAD = 0;

    // time loop
    for (int iTime = 0; iTime < nTime; ++iTime){
        MPI_Send(&currGrid[1][0], nColsWithGhost, MPI_INT, upperNeighbor, 0,
                MPI_COMM_WORLD);
        MPI_Send(&currGrid[nRowsLocal][0], nColsWithGhost, MPI_INT,
                lowerNeighbor, 0, MPI_COMM_WORLD);

        MPI_Recv(&currGrid[nRowsLocal + 1][0], nColsWithGhost, MPI_INT,
                lowerNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&currGrid[0][0], nColsWithGhost, MPI_INT, upperNeighbor, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int iRow = 0; iRow < nRowsLocalWithGhost; ++iRow){
            currGrid[iRow][0] = currGrid[iRow][nCols];
            currGrid[iRow][nCols + 1] = currGrid[iRow][1];
        }

        // display
        if (mpirank != mpiroot) {
            for (int iRow = 1; iRow <= nRowsLocal; ++iRow){
                MPI_Send(&currGrid[iRow][1], nCols, MPI_INT, mpiroot, 0,
                        MPI_COMM_WORLD);
            }
        }

        if (mpirank == mpiroot) {
            cout << "iTime: " << iTime << "-----------------------" << endl;

            // print its own data
            for (int iRow = 1; iRow <= nRowsLocal; ++iRow){
                for (int iCol = 1; iCol <= nCols; ++iCol){
                    cout << currGrid[iRow][iCol] << " ";
                }
                cout << endl;
            }

            // receive from other ranks
            for (int sourceRank = 1; sourceRank < mpisize; ++sourceRank){
                int nRecv = nRows / mpisize;
                // last rank
                if (sourceRank == mpisize - 1){
                    nRecv += nRows % mpisize;
                }

                // FIX: nRecv meaning
                vector<int> buff(nCols, 0);
                for (int iRecv = 0; iRecv < nRecv; ++iRecv){
                    MPI_Recv(&buff[0], nCols, MPI_INT, sourceRank, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int i : buff){
                        cout << i << " ";
                    }
                    cout << endl;
                }
            }
        }

        // update the grid
        for (int iRow = 1; iRow <= nRowsLocal; ++iRow) {
            for (int iCol = 1; iCol <= nCols; ++iCol){
                int nAliveNeighbors = 0;

                for (int jRow = iRow - 1; jRow <= iRow + 1; ++jRow){
                    for (int jCol = iCol - 1; jCol <= iCol + 1; ++jCol){
                        if ( (jRow != iRow || jCol != iCol) &&
                        currGrid[jRow][jCol] == ALIVE){
                            ++nAliveNeighbors;
                        }
                    }
                }

                if (nAliveNeighbors < 2) {
                    nextGrid[iRow][iCol] = DEAD;
                }

                if (currGrid[iRow][iCol] == ALIVE && (nAliveNeighbors == 2 ||
                nAliveNeighbors == 3) ) {
                    nextGrid[iRow][iCol] = ALIVE;
                }

                if (nAliveNeighbors > 3){
                    nextGrid[iRow][iCol] = DEAD;
                }

                if (currGrid[iRow][iCol] == DEAD && nAliveNeighbors == 3){
                    nextGrid[iRow][iCol] = ALIVE;
                }
            }
        }

        // swap the grids
//        currGrid.swap(nextGrid);
        for (int iRow = 1; iRow <= nRowsLocal; ++iRow){
            for (int iCol = 1; iCol <= nCols; ++ iCol){
                currGrid[iRow][iCol] = nextGrid[iRow][iCol];
            }
        }

    }

    MPI_Comm old_comm, new_comm;
    int ndims, reorder, periods[2], dim_size[2];
    MPI_Finalize();
    return 0;
}