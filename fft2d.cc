// Distributed two-dimensional Discrete FFT transform
// YOUR NAME HERE: Mayank Parasar
// ECE8893 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <signal.h>
#include <math.h> // This is included because W = cos(2*M_PI/width) - j*sin(2*M_PI/width)
#include <mpi.h>
#include <assert.h>

#include "Complex.h"
#include "InputImage.h"
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
using namespace std;

void Transform1D(Complex* h, int w, Complex* H, int row_per_cpu);
void invTransform1D(Complex* h, int w, Complex* H, int row_per_cpu);

void Transform2D(const char* inputFN) //<--->Transform2D(fn.c_str())
{ 

  InputImage image(inputFN);  // Create the helper object for reading the image [inputFN = Tower.txt; 
			      // Constructor contains the pointer to file to read from it and populate 
			      // the Complex structure for the input file.]
  // Step 1:
  // Use the InputImage object to read in the Tower.txt file and
  // find the width/height of the input image.
  // Your code here, steps 2-9
	// Author: Mayank Parasar
	// gtid: 902978052
	// contact: mparasar3@gatech.edu
	int numtask, rank; 
	int width = image.GetWidth();
	int height = image.GetHeight();
	 // Step 2:
	 // Use MPI to find how many CPUs in total, and which one this 
	 // CPU
	MPI_Comm_size(MPI_COMM_WORLD, &numtask);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	 // Step 3: 
	 // Allocate an array of Complex object of sufficient size to
	 // hold the 2d DFT results (size is width * height)
	assert((height % numtask) == 0); 
	 // Always allocate memory to a pointer
	Complex *out_data = new Complex[width*height]; 
				    // Pointer to an object of class type
					// Complex, used to store 1D and 
					// eventually 2D transform.
					// Complex is basically a class describing
					// One complex number. So for each complex
					// number we would need a separate object
					// hence is the array of object of class 'Complex'	
	// Step 4: Obtain a pointer to the Complex 1d array of input data.
	// Always allocate memory to a pointer
	// 'inp_data' is Pointer to an object of class 'Complex'
	Complex *inp_data = new Complex[width*height]; 
	//class Complex *out_data;
  	inp_data = image.GetImageData(); 
  					// GetImageData() returns a one dimentional array
					// of type 'Complex' representing the original
					// time-domain input values 
	// Step 5: Do the individual 1D transforms on the rows assigned to your CPUv
	int row_per_cpu, start_row;
	row_per_cpu = height/numtask; // Row count to each CPU 
	start_row = rank * row_per_cpu; // Starting row for each CPU 
	
	// Pass array, with proper indecies to a function, rather than pointer 	
	
	// Each CPU will perform it's own 1D-transform
	Transform1D(&inp_data[rank * row_per_cpu * width], width, out_data, row_per_cpu);

	int rc;
	MPI_Status status;
	MPI_Request pReq;
	if (rank == 0) {
	// Rank 0 should recieve from all other ranks' 1D
	// If it's CPU 0 then it should receive from all other CPUs
	// Each receive should then be as many times as there are
	// other CPUs
	for (int i = 1; i < numtask; i++) {
		rc =  MPI_Irecv((&out_data[i*width*row_per_cpu]),	// buffer(pointer to where data is to be put in buffer?)
				row_per_cpu * width * sizeof(Complex),	// Count
				MPI_CHAR,				// type
				i,					// source
				0,		 			// tag
				MPI_COMM_WORLD,				// comm
				//&status
				&pReq	// request
				);
		if (rc != MPI_SUCCESS)  {
              		cout << "Rank " << rank
               	    	<< " send failed, rc " << rc << endl;
              		MPI_Finalize();
              		exit(1);
		}
		MPI_Wait(&pReq,&status);
		}

	// Write this data to a file
	string fo("MyAfter1d.txt");
	image.SaveImageData(fo.c_str(), out_data,
                		               256, 256);
       } else { 		
// 	All other CPUs will send to CPU 0
//	Each CPU will send it's 'out_data', thus each CPU needs to call
//	this API only once!
		rc = MPI_Send(out_data,			// Buffer
			row_per_cpu * width * sizeof(Complex),	// Count
			MPI_CHAR,				// Datatype
			0,					// Dest
			0,					// Tag
			MPI_COMM_WORLD				// Comm
			);

	if (rc != MPI_SUCCESS)  {
            	cout << "Rank " << rank
               	<< " send failed, rc " << rc << endl;
	       	MPI_Finalize();
              	exit(1);
	}
}

	/***************2-D Transform******************/

	Complex *tran_inp_data = new Complex[width*height]; 
	Complex *tran_out_data = new Complex[width*height]; 
	Complex *final_out_data = new Complex[width*height]; 
	if (rank == 0) {
// Transpose: Rank 0 will do the transpose
		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++)
				tran_inp_data[c*width + r] = out_data[r*width + c];
// Send to all CPUs all of the data: for_loop
		for(int i = 1; i < numtask; i++) {
			rc = MPI_Send(tran_inp_data,			// Buffer
				height * width * sizeof(Complex),	// Count
				MPI_CHAR,				// Datatype
				i,					// Dest
				0,					// Tag
				MPI_COMM_WORLD				// Comm
				);
			if (rc != MPI_SUCCESS)  {
            			cout << "Rank " << rank
               			<< " send failed, rc " << rc << endl;
	       			MPI_Finalize();
              			exit(1);
			}
		}
	} else {
// All CPUs will receive
		rc =  MPI_Recv(tran_inp_data,	// buffer(pointer to where data is to be put in buffer?)
			height * width * sizeof(Complex),	// Count
			MPI_CHAR,				// type
			0,					// source
			0,		 			// tag
			MPI_COMM_WORLD,				// comm
			&status
			//&pReq	// request
			);
		if (rc != MPI_SUCCESS)  {
              		cout << "Rank " << rank
               	    	<< " send failed, rc " << rc << endl;
              		MPI_Finalize();
              		exit(1);
			}
	}

// All CPUs, including CPU 0, will call 'Transform1D' again

	Transform1D(&tran_inp_data[rank * row_per_cpu * width], width, tran_out_data, row_per_cpu);

	if (rank == 0) {
// Rank 0 receive from all CPUs, again
	for (int i = 1; i < numtask; i++) {
		rc =  MPI_Recv((&tran_out_data[i*width*row_per_cpu]),	// buffer(pointer to where data is to be put in buffer?)
				row_per_cpu * width * sizeof(Complex),	// Count
				MPI_CHAR,				// type
				i,					// source
				0,		 			// tag
				MPI_COMM_WORLD,				// comm
				&status
				);
		if (rc != MPI_SUCCESS)  {
              		cout << "Rank " << rank
               	    	<< " send failed, rc " << rc << endl;
              		MPI_Finalize();
              		exit(1);
			}
		}
// To the final transpose of the tran_out_data->final_out_data
		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++)
				final_out_data[c*width + r] = tran_out_data[r*width + c];

// Write this data to a file to compare after2d.txt
	string f("MyAfter2d.txt");
	image.SaveImageData(f.c_str(), final_out_data,
                		               256, 256);
	} else {
// All other CPUs will send the data
		rc = MPI_Send(tran_out_data,			// Buffer
			row_per_cpu * width * sizeof(Complex),	// Count
			MPI_CHAR,				// Datatype
			0,					// Dest
			0,					// Tag
			MPI_COMM_WORLD				// Comm
			);
		if (rc != MPI_SUCCESS)  {
              		cout << "Rank " << rank
               	    	<< " send failed, rc " << rc << endl;
              		MPI_Finalize();
              		exit(1);
		}
	}


		/**********InverseDFT**************/
	Complex *inv_inp_data = new Complex[width*height];
  if (rank == 0)
	inv_inp_data = final_out_data; //'final_out_data' is only with CPU0, contains the 2-D DFT 
	Complex *inv_out_data = new Complex[width*height]; 
	Complex *tran_inv_out_data = new Complex[width*height]; 
	Complex *tran_final_inv_out_data = new Complex[width*height];
	Complex *final_inv_out_data = new Complex[width*height];
 
// CPU0 Send the whole file to all CPUs
	if (rank == 0) {
		for(int i = 1; i < numtask; i++) {
			rc = MPI_Send(inv_inp_data,			// Buffer
				height * width * sizeof(Complex),	// Count
				MPI_CHAR,				// Datatype
				i,					// Dest
				0,					// Tag
				MPI_COMM_WORLD				// Comm
				);
		if (rc != MPI_SUCCESS)  {
   			cout << "Rank " << rank
       			<< " send failed, rc " << rc << endl;
      			MPI_Finalize();
       			exit(1);
		}
		}
	} else {
// All CPUs will receive the whole file
		rc =  MPI_Recv(inv_inp_data,	// buffer(pointer to where data is to be put in buffer?)
			height * width * sizeof(Complex),	// Count
			MPI_CHAR,				// type
			0,					// source
			0,		 			// tag
			MPI_COMM_WORLD,				// comm
			&status
			);
		if (rc != MPI_SUCCESS)  {
        	cout << "Rank " << rank
           	<< " send failed, rc " << rc << endl;
           	MPI_Finalize();
           	exit(1);
		}
	}

// All CPUs, including CPU0 will call invTransform1D
// with their allotted size of received file

	invTransform1D(&inv_inp_data[rank * row_per_cpu * width], width, inv_out_data, row_per_cpu);

// CPU0 will receive from all CPUs
if (rank == 0) {
	for (int i = 1; i < numtask; i++) {
		rc =  MPI_Recv((&inv_out_data[i*width*row_per_cpu]),	// buffer(pointer to where data is to be put in buffer?)
				row_per_cpu * width * sizeof(Complex),	// Count
				MPI_CHAR,				// type
				i,					// source
				0,		 			// tag
				MPI_COMM_WORLD,				// comm
				&status
				);
		if (rc != MPI_SUCCESS)  {
           cout << "Rank " << rank
           << " send failed, rc " << rc << endl;
           MPI_Finalize();
           exit(1);
			}
		}
// the transpose of the inv_out_data->tran_inv_out_data
		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++)
				tran_inv_out_data[c*width + r] = inv_out_data[r*width + c];
	} else {
// All other CPUs will send the data
	rc = MPI_Send(inv_out_data,			// Buffer
		row_per_cpu * width * sizeof(Complex),	// Count
		MPI_CHAR,				// Datatype
		0,					// Dest
		0,					// Tag
		MPI_COMM_WORLD				// Comm
		);
	if (rc != MPI_SUCCESS)  {
		cout << "Rank " << rank
        	<< " send failed, rc " << rc << endl;
    	MPI_Finalize();
    	exit(1);
		}
	}
// CPU 0 again send the whole file to all the CPUs
if (rank == 0) {
	for(int i = 1; i < numtask; i++) {
		rc = MPI_Send(tran_inv_out_data,			// Buffer
			height * width * sizeof(Complex),	// Count
			MPI_CHAR,				// Datatype
			i,					// Dest
			0,					// Tag
			MPI_COMM_WORLD				// Comm
			);
	if (rc != MPI_SUCCESS)  {
		cout << "Rank " << rank
  			<< " send failed, rc " << rc << endl;
   			MPI_Finalize();
   			exit(1);
		}
	}
} else {
// All CPUs will receive the whole file
		rc =  MPI_Recv(tran_inv_out_data,	// buffer(pointer to where data is to be put in buffer?)
			height * width * sizeof(Complex),	// Count
			MPI_CHAR,				// type
			0,					// source
			0,		 			// tag
			MPI_COMM_WORLD,				// comm
			&status
			);
		if (rc != MPI_SUCCESS)  {
        	cout << "Rank " << rank
           	<< " send failed, rc " << rc << endl;
           	MPI_Finalize();
           	exit(1);
		}
	}
// All CPUs, including CPU0 will call Transform1D
// with their allotted size of received file

	invTransform1D(&tran_inv_out_data[rank * row_per_cpu * width], width, tran_final_inv_out_data, row_per_cpu);

// CPU0 will receive from all CPUs
if (rank == 0) {
	for (int i = 1; i < numtask; i++) {
		rc =  MPI_Recv((&tran_final_inv_out_data[i*width*row_per_cpu]),	// buffer(pointer to where data is to be put in buffer?)
				row_per_cpu * width * sizeof(Complex),	// Count
				MPI_CHAR,				// type
				i,					// source
				0,		 			// tag
				MPI_COMM_WORLD,				// comm
				&status
				);
		if (rc != MPI_SUCCESS)  {
           cout << "Rank " << rank
           << " send failed, rc " << rc << endl;
           MPI_Finalize();
           exit(1);
			}
		}
		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++)
				final_inv_out_data[c*width + r] = tran_final_inv_out_data[r*width + c];
// CPU0 writes the result in a new file... GameOver!!!
// Write this data to a file to compare after2d.txt
	string finv("MyAfterInverse.txt");
	image.SaveImageData(finv.c_str(), final_inv_out_data,
                		               256, 256);
	} else {
// All other CPUs will send the data
	rc = MPI_Send(tran_final_inv_out_data,			// Buffer
		row_per_cpu * width * sizeof(Complex),	// Count
		MPI_CHAR,				// Datatype
		0,					// Dest
		0,					// Tag
		MPI_COMM_WORLD				// Comm
		);
	if (rc != MPI_SUCCESS)  {
		cout << "Rank " << rank
        	<< " send failed, rc " << rc << endl;
    	MPI_Finalize();
    	exit(1);
		}
	}
}

void Transform1D(Complex* h, int w, Complex* H, int row_per_cpu)
{
	  // Implement a simple 1-d DFT using the double summation equation
	  // given in the assignment handout.  h is the time-domain input
	  // data, w is the width (N), and H is the output array.

	for (int k = 0; k < row_per_cpu; k++) {
	  for (int i = 0; i < w; i++) {	// Each element of H
  		Complex W;
		H[(i + k*w)] = Complex(0,0);
		 for (int j = 0; j < w; j++) {	// do multiply and add on all other elements
			W.real = cos((2 * M_PI * i * j)/w); W.imag = -sin((2 * M_PI * i * j)/w);
			H[(i + k*w)] = H[(i + k*w)] + (W * h[(j + k*w)]);
			}
		if (fabs(H[i].imag) < 1e-10) H[i].imag = 0;
		if (fabs(H[i].real) < 1e-10) H[i].real = 0;
		}
	}
}

void invTransform1D(Complex* h, int w, Complex* H, int row_per_cpu)
{
	  // Implement a simple 1-d DFT using the double summation equation
	  // given in the assignment handout.  h is the time-domain input
	  // data, w is the width (N), and H is the output array.

	for (int k = 0; k < row_per_cpu; k++) {
	  for (int i = 0; i < w; i++) {
  		Complex W;
		H[(i + k*w)] = Complex(0,0);
		 for (int j = 0; j < w; j++) {	
			W.real = cos((2 * M_PI * i * j)/w); W.imag = sin((2 * M_PI * i * j)/w);
			H[(i + k*w)] = H[(i + k*w)] + (W * h[(j + k*w)]);
			}
		// Divide by width(N) in the function itself
		H[(i + k*w)].real = H[(i + k*w)].real/w; H[(i + k*w)].imag = H[(i + k*w)].imag/w;
		if (fabs(H[i].imag) < 1e-10) H[i].imag = 0;
		if (fabs(H[i].real) < 1e-10) H[i].real = 0;
		}
	}
}


int main(int argc, char** argv)
{
  int rc;

  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS ) {
	printf("Error starting MPI program. Termination. \n");
	MPI_Abort(MPI_COMM_WORLD, rc);
	}
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.

  MPI_Finalize();

  return 0;
}  
  

  
