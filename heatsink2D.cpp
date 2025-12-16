/*
ASTE-404 HW6 initial file
*/

#include <iostream>         // for screen output
#include <fstream>          // for file writing
#include <cmath>

// function prototypes
bool solveGS_SOR(double *a, double *b, double *c, double *d, double *e, double *g, double *T, int ni, int nj);
double NomadContour(double x);

int main() {

    // Mesh Sizing for Chamber Profile
    int ni = 960;   // number of nodes
    int nj = 213;

    double x0 = 0.0;  // origin
    double y0 = 0.0;

    double dx = 0.01 * 0.0254;  // cell spacing converted from in to m
    double dy = 0.01 * 0.0254;

    int nn = ni*nj;  // total number of nodes

    // Define solid vs. gas nodes
    char *solid = new char[nn];     // initial array for 0 = gas, 1 = solid
    for (int n = 0; n < nn; n++){   // clear data
        solid[n] = 0;
    }

    for (int j=0; j<nj; j++){
        for (int i=0; i<ni; i++){
            double y = y0 + j * dy;
            double x = x0 + i * dx;
            double h = 0;
                        
            

        }
    }


    // Memory Allocation
    double *T = new double[nn];   // temperature array
    for (int n=0; n<nn; n++){     // clear data
        T[n] = 0;
    }

    double *a = new double[nn];   // matrix coefficients
    double *b = new double[nn];
    double *c = new double[nn];
    double *d = new double[nn];
    double *e = new double[nn];
    double *g = new double[nn];   // RHS 
    for (int n=0; n<nn; n++){     // clear data
        a[n] = 0;
        b[n] = 0;
        c[n] = 0;
        d[n] = 0;
        e[n] = 0;
        g[n] = 0;
    }
    
    // Set Matrix Values
    double dxsqr = dx*dx;
    double dysqr = dy*dy;
    for (int j=0; j<nj; j++){
        for (int i=0; i<ni; i++){
            // initialize unknown indexing variable
            int n = j*ni + i;

            // start with internal nodes and then overwrite boundary nodes
            a[n] = e[n] = 1/(dysqr); // assign multiple variables
            b[n] = d[n] = 1/(dxsqr);
            c[n] = -2/(dxsqr)-2/(dysqr);

            // boundaries
            if (i==0){
                // zero Neumann on x min
                b[n] = 0;
                d[n] = 2/(dxsqr);
            }
            else if (i==ni-1){
                // zero Neumann on x max
                b[n] = 2/(dxsqr);
                d[n] = 0;
            }
            if (j==0){
                // Dirichlet on y min
                a[n] = b[n] = d[n] = e[n] = 0;
                c[n] = 1;
                g[n] = 300;
                continue; // node fully determined now so skip to next
            }
            else if (j==nj-1){
                // zero neumann on y max
                a[n] = 2/(dysqr);
                e[n] = 0;
            }

        }

    } // end matrix loops

    // set random internal Dirichlet points
    for (int s = 0; s<20; s++){
        int i = (int) (1 + rnd()*(ni-2)); // integer in [1,ni-2]
        int j = (int) (1 + rnd()*(nj-2)); // integer in [1,nj-2]
        int n = j*ni + i;

        //clear all row data
        a[n] = b[n] = c[n] = d[n] = e[n] = g[n] = 0;

        // make Dirichlet node
        c[n] = 1.0; // main diagonal
        g[n] = rnd()*500 + 300; // random value in [300,800)
    }

    //call solver
    solveGS_SOR(a,b,c,d,e,g,T,ni,nj);

    // free memory
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] e;
    delete[] g;

    // output to VTI file
    std::ofstream out("field.vti");

    out<<"<VTKFile type=\"ImageData\">\n";
    out<<"<ImageData WholeExtent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<0<<"\"";
    out<<" Origin=\""<<x0<<" "<<y0<<" "<<0.0<<"\"\n";
    out<<" Spacing=\""<<dx<<" " <<dy<<" "<<0.0<<"\">\n";
    out<<"<Piece Extent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<0<<"\">\n";
    out<<"<PointData>\n";

    out<<"<DataArray Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
    for (int n=0;n<nn;n++) out<<T[n]<<" ";
    out<<"\n</DataArray>\n";

    out<<"</PointData>\n";
    out<<"</Piece>\n";
    out<<"</ImageData>\n";
    out<<"</VTKFile>\n";

    // free remaining temperature array
    delete[] T;
    return 0;	// normal exit
}

bool solveGS_SOR(double *a, double *b, double *c, double *d, double *e, double *g, double *T, int ni, int nj){

    int nn = ni*nj;

    // Solve Matrix System
    const double w = 1.4; // SOR relaxation factor
    for (int it=0; it<10000; it++){ // solver iteration
        for (int n=0; n<nn; n++){ // loop over all nodes
            // compute row n * vec(T) product but skip over diagonal term
            // we only include non-zero entries to avoid trying to access out of bound T entries
            double sum = 0;
            if (a[n]!=0) sum += a[n]*T[n-ni]; // T[i,j-1] term
            if (b[n]!=0) sum += b[n]*T[n-1]; // T[i-1,j] term
            if (d[n]!=0) sum += d[n]*T[n+1]; // T[i+1,j] term
            if (e[n]!=0) sum += e[n]*T[n+ni]; // T[i,j+1] term
            double T_star = (g[n] - sum) / c[n]; // new estimate for T[i,j]

            // perform SOR step here to update T[n]
            T[n] = T[n] + w*(T_star - T[n]);
        }

        // Check for Convergence
        if (it%50==0){
            double r2_sum = 0;
            for (int n=0; n<nn; n++){
                double sum = 0;
                if (a[n]!=0) sum += a[n]*T[n-ni]; // T[i,j-1] term
                if (b[n]!=0) sum += b[n]*T[n-1]; // T[i-1,j] term
                sum += c[n]*T[n];
                if (d[n]!=0) sum += d[n]*T[n+1]; // T[i+1,j] term
                if (e[n]!=0) sum += e[n]*T[n+ni]; // T[i,j+1] term

                double r = g[n] - sum;
                r2_sum += r*r;
            }

            // compute avg error
            double L2 = sqrt(r2_sum/nn);

            std::cout<<"solve iteration: "<<it<<", L2 norm: "<<L2<<std::endl;
            if (L2<1e-6) return true; // break out of loop when converged
        }
    }
    return false;
}

double NomadContour(double x){
    x = x / 0.0254; // convert to inches
     
    double h = 0.0;
    double num = 0.0;

    if (x < 5.94) {
        h = 0.82 * 0.0254;
    }
    else if (x < 6.56) {
        num = 0.93 - (x - 5.94)^2
        h = -0.14 + sqrt(0.93 - (x - 5.94)^2);

    }
}
