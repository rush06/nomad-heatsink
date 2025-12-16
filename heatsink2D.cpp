/*
ASTE-404 HW6 initial file
*/

#include <iostream>         // for screen output
#include <fstream>          // for file writing
#include <cmath>
#include <vector>

// function prototypes
bool solveGS_SOR(double *a, double *b, double *c, double *d, double *e, double *g, double *T, int ni, int nj, const char* solid);
double NomadContour(double x);
double spline(double x);

int main() {

    // Mesh Sizing for Chamber Profile
    int ni = 960;   // number of nodes
    int nj = 213;

    double x0 = 0.0;  // origin
    double y0 = 0.0;

    double dx = 0.01;  // cell spacing converted from in to m
    double dy = 0.01;

    int nn = ni*nj;  // total number of nodes

    // Memory Allocation
    double *T = new double[nn];   // temperature array
    for (int n=0; n<nn; n++){     // clear data and set initial guess
        T[n] = 0;
    }

    // Define solid vs. gas nodes
    char *solid = new char[nn];     // initial array for 0 = gas, 1 = solid
    for (int n = 0; n < nn; n++){   // clear data
        solid[n] = 0;
    }

    for (int j=0; j<nj; j++){       // loop to set solid and gas nodes
        for (int i=0; i<ni; i++){
            int n = j * ni + i;
            double y = y0 + j * dy;
            double x = x0 + i * dx;
            double h = NomadContour(x);
            
            if (y < h){
                solid[n] = 0;       // gas nodes
                T[n] = 3670.0;
            }
            else {
                solid[n] = 1;       // solid nodes
                T[n] = 300.0;
            }

            // Boundaries on top are ambient gas
            if (j == nj-1){solid[n] = 0;}
        }
    }

    double *a = new double[nn];   // matrix coefficients
    double *b = new double[nn];
    double *c = new double[nn];
    double *d = new double[nn];
    double *e = new double[nn];
    double *g = new double[nn];   // RHS 
    double *jWall = new double[nn]; // used for finding wall indices
    for (int n=0; n<nn; n++){     // clear data
        a[n] = 0;
        b[n] = 0;
        c[n] = 0;
        d[n] = 0;
        e[n] = 0;
        g[n] = 0;
        jWall[n] = 0;
    }
    
    // get hot wall indices
    for (int i = 0; i < ni; ++i) {
        double xP = x0 + i * dx;
        double h = NomadContour(xP);        // inner wall y
        
        int jw = int((h - y0) / dy + 0.5);  // closest grid index to the wall
        if (jw < 0)      jw = 0;
        if (jw > nj - 1) jw = nj - 1;

        jWall[i] = jw;
    }

    // Set Matrix Values
    double k = 401.0;       // heat transfer coefficient for copper
    double dxsqr = dx*dx;
    double dysqr = dy*dy;
    for (int j=0; j<nj; j++){
        for (int i=0; i<ni; i++){
            // initialize unknown indexing variable
            int n = j*ni + i;

            if (solid[n] == 0){
                a[n] = b[n] = d[n] = e[n] = 0.0;
                c[n] = 1.0;
                continue;
            }

            // start with internal nodes and then overwrite boundary nodes
            a[n] = e[n] = k/(dysqr); // assign multiple variables
            b[n] = d[n] = k/(dxsqr);
            c[n] = -2*k/(dxsqr) -2*k/(dysqr);

            // Convection BCs
            double xP = x0 + i * dx;
            double yP = y0 + j * dy;
            if (j == jWall[i] || (xP < 5.94 && yP == 0.82)){     // lower BC
                double yin = NomadContour(xP);     // inner wall y at this x
                double delta = yP - yin;           // distance from node to wall in -y direction

                // avoid grid mismatch
                if (delta < 1e-12) delta = 1e-12;
                if (delta > dy)    delta = dy;
                double f = delta / dy;

                // beta factor
                double coef = k / (dy * dy);
                double h_hot = 15000;           // assumption - can be added later with Bartz Equation
                double T_gas = 3670;            // assumption - can be added later
                double beta = h_hot / (k/dy + h_hot * f);

                // modify neighbor nodes
                a[n] = 0.0;
                c[n] += k / (dy * dy);
                c[n] -= coef * beta;
                g[n] -= coef * beta * T_gas;
            }
            if (j == nj-2){    // upper BC
                double yout = 2.12;                 // outer wall y
                double delta = yout - yP;           // distance from node to wall in +y direction

                // avoid grid mismatch
                if (delta < 1e-12) delta = 1e-12;
                if (delta > dy)    delta = dy;
                double f = delta / dy;

                // beta factor
                double coef = k / (dy * dy);
                double h_amb = 20;       // assumption for ambient air
                double T_amb = 300;        // assumption for ambient air
                double beta = h_amb / (k/dy + h_amb * f);

                // modify neighbor nodes
                e[n] = 0.0;
                c[n] += k / (dy * dy);
                c[n] -= coef * beta;
                g[n] -= coef * beta * T_amb;
            }

            // Neumann BCs
            if (i==0){
                // zero Neumann on x min
                b[n] = 0;
                d[n] = 2*k/(dxsqr);
            }
            else if (i==ni-1){
                // zero Neumann on x max
                b[n] = 2*k/(dxsqr);
                d[n] = 0;
            }
        }

    } // end matrix loops

    //call solver
    solveGS_SOR(a,b,c,d,e,g,T,ni,nj,solid);

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

bool solveGS_SOR(double *a, double *b, double *c, double *d, double *e, double *g, double *T, int ni, int nj, const char* solid){

    int nn = ni*nj;
    const double w = 1.0; // SOR relaxation factor

    // Solve Matrix System
    for (int it=0; it<10000; it++){ // solver iteration
        for (int n=0; n<nn; n++){ // loop over all nodes
            // skip gas nodes
            if (!solid[n]){continue;}

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

        // convergence check
        if (it%50==0){
            double r2_sum = 0;
            int skipped = 0;
            for (int n=0; n<nn; n++){
                // skip gas nodes
                if (!solid[n]){skipped++;continue;}

                double sum = 0;
                if (a[n]!=0) sum += a[n]*T[n-ni]; // T[i,j-1] term
                if (b[n]!=0) sum += b[n]*T[n-1]; // T[i-1,j] term
                sum += c[n]*T[n];
                if (d[n]!=0) sum += d[n]*T[n+1]; // T[i+1,j] term
                if (e[n]!=0) sum += e[n]*T[n+ni]; // T[i,j+1] term

                double r = g[n] - sum;
                r2_sum += r*r;
            }

            // compute avg error ONLY for solid nodes
            double L2 = sqrt(r2_sum/(nn-skipped));

            std::cout<<"solve iteration: "<<it<<", L2 norm: "<<L2<<std::endl;
            if (L2<1e-6) return true; // break out of loop when converged
        }
    }
    return false;
}

double NomadContour(double x){
    double h = 0.0;

    if (x < 5.94) {
        h = 0.82;
    }
    else if (x < 6.56) {
        h = -0.14 + sqrt(0.93 - (x - 5.94)*(x - 5.94));
    }
    else if (x < 7.00){
        h = 0.6 - 0.84 * (x - 6.56);
    }
    else if (x < 7.62){
        h = 0.9637 - sqrt(0.9637*0.9637 - (x-7.62)*(x-7.62));
    }
    else if (x < 7.71){
        h = 0.245 - sqrt(0.245*0.245 - (x-7.62)*(x-7.62));
    }
    else if (x <= 9.6){
        double offset = 0.017129 - spline(7.71);
        h = spline(x) + offset;
    }

    return h;
}

double spline(double x){
    x = x - 7.71;
    double h = 0.0;

    // lookup table from engine spline
    std::vector<double> xs = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.89};
    std::vector<double> hs = {0, 0.061, 0.103, 0.145, 0.187, 0.228, 0.269, 0.309, 0.349, 0.387, 0.424, 0.459, 0.493, 0.525, 0.555, 0.583, 0.608, 0.630, 0.650, 0.664}; 

    // interpolate between points
    for (int i=1; i<xs.size(); i++){
        if (x < xs[i]){
            h = (hs[i] - hs[i-1])/(xs[i] - xs[i-1]) * (x - xs[i-1]) + hs[i-1];
            break;
        }
    }

    return h;
}
