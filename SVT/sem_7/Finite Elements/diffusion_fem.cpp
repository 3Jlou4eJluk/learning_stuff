#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

// // Corresponds to tensor
// // [ 1  0 ]
// // [ 0 10 ]
// // rotated by M_PI/6
// const double Dxx = 3.25;
// const double Dyy = -0.433013;
// const double Dxy = 0.25;

const double dx = 1.0;
const double dy = 1.0;
const double dxy = 0.0;
const double pi = 3.1415926535898;
const double a = 1;

double nodeDist(const Node &n1, const Node &n2)
{
    double x1[3], x2[3];
    n1.Centroid(x1);
    n2.Centroid(x2);
    return sqrt(
            (  x1[0] - x2[0])*(x1[0] - x2[0])
            + (x1[1] - x2[1])*(x1[1] - x2[1])
            + (x1[2] - x2[2])*(x1[2] - x2[2])
    );
}

double cellDiam(const Cell &c)
{
    ElementArray<Node> nodes = c.getNodes();
    unsigned m = static_cast<unsigned>(nodes.size());
    double diam = 0.;
    for(unsigned i = 0; i < m; i++){
        for(unsigned j = 0; j < m; j++){
            diam = max(diam, nodeDist(nodes[i], nodes[j]));
        }
    }
    return diam;
}

void main_mesh_diam(Mesh &m)
{
    double diam = 0.0;
    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        diam = max(diam, cellDiam(icell->self()));
    }
    cout << "Mesh diameter is " << diam << endl;
}

double C(double x, double y)
{
    //return 1;
    //return x;
	return sin(a*x) * sin(a*y);
}

double source(double x, double y)
{
	//return 0;
	return -a*a * (2.*dxy * cos(a*x)*cos(a*y) - (dx+dy) * sin(a*x)*sin(a*y));
}

enum BoundCondType
{
	BC_DIR = 1,
	BC_NEUM = 2
};

// Class including everything needed
class Problem
{
private:
	/// Mesh
	Mesh &m;
	// =========== Tags =============
	/// Solution tag: 1 real value per node
	Tag tagConc;
	/// Diffusion tensor tag: 3 real values (Dx, Dy, Dxy) per cell
	Tag tagD;
	/// Boundary condition type tag: 1 integer value per node
	Tag tagBCtype;
	/// Boundary condition value tag: 1 real value per node, sparse on nodes
	Tag tagBCval;
	/// Right-hand side tag: 1 real value per node, sparse on nodes
	Tag tagSource;
	/// Analytical solution tag: 1 real value per node
	Tag tagConcAn;
	/// Global index tag: 1 integer value per node
	Tag tagGlobInd;
    /// Transmissibility: 1 real value per edge
    Tag tagTrans;

    Tag tagFlux;

	// =========== Tag names ===========
	const string tagNameConc = "Concentration";
	const string tagNameD = "Diffusion_tensor";
	const string tagNameBCtype = "BC_type";
	const string tagNameBCval = "BC_value";
	const string tagNameSource = "Source";
	const string tagNameConcAn = "Concentration_analytical";
	const string tagNameGlobInd = "Global_Index";
    const string tagNameTrans = "Transmissibility";
    const string tagNameFlux = "Flux";

	// =========== Markers
	/// Marker for Dirichlet nodes
	MarkerType mrkDirNode;
	/// Number of Dirichlet nodes
	unsigned numDirNodes;
public:
	Problem(Mesh &m_);
	~Problem();
	void initProblem();
	void assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs);
	void assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc);
	void run();
    double C_norm();
    double L2_norm();
};

Problem::Problem(Mesh &m_) : m(m_)
{
}

Problem::~Problem()
{

}

void Problem::initProblem()
{
	// Init tags
	tagConc = m.CreateTag(tagNameConc, DATA_REAL, CELL, NONE, 1);
	tagD = m.CreateTag(tagNameD, DATA_REAL, CELL, NONE, 3);
	tagBCtype = m.CreateTag(tagNameBCtype, DATA_INTEGER, FACE, FACE, 1);
	tagBCval = m.CreateTag(tagNameBCval, DATA_REAL, FACE, FACE, 1);
	tagSource =  m.CreateTag(tagNameSource, DATA_REAL, CELL, NONE, 1);
	tagConcAn =  m.CreateTag(tagNameConcAn, DATA_REAL, CELL, NONE, 1);
	tagGlobInd = m.CreateTag(tagNameGlobInd, DATA_INTEGER, CELL, NONE, 1);
    tagTrans = m.CreateTag(tagNameTrans, DATA_REAL, FACE, NONE, 1);
    tagFlux = m.CreateTag(tagNameFlux, DATA_REAL, CELL, NONE, 1);

	// Cell
    int glob_ind = 0;
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();

        double xn[3];
        c.Barycenter(xn);
		
		c.RealArray(tagD)[0] = dx;
		c.RealArray(tagD)[1] = dy;
		c.RealArray(tagD)[2] = dxy;

        c.Real(tagConcAn) =  C(xn[0], xn[1]);
        c.Real(tagSource) = source(xn[0], xn[1]);

        c.Integer(tagGlobInd) = glob_ind;
        glob_ind++;
	}

	// Face loop
	mrkDirNode = m.CreateMarker();
	numDirNodes = 0;
	for(Mesh::iteratorFace face = m.BeginFace(); face != m.EndFace(); face++){
		Face f = face->getAsFace();

		if(f.Boundary()){
			f.SetMarker(mrkDirNode);
			f.Integer(tagBCtype) = BC_DIR;
			numDirNodes++;

            ElementArray<Cell> c = f->getCells();
			f.Real(tagBCval) = c[0].Real(tagConcAn);
		}
	}
}

double basis_func(const Cell &c, const Node &n, double x_, double y_)
{
    ElementArray<Node> nodes = c.getNodes();
	unsigned n_ind = 0;
	double x[3];
	double y[3];
	for(unsigned i = 0; i < 3; i++){
		if(n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}
	
	if(n_ind == 0){
		return ((x_   - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y_   - y[2])) /
			   ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2]));
	}
	else if(n_ind == 1){
		return ((x_   - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y_   - y[2])) /
			   ((x[1] - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y[1] - y[2]));
	}
	else if(n_ind == 2){
		return ((x_   - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y_   - y[0])) /
			   ((x[2] - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y[2] - y[0]));
	}
	else{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
}

rMatrix basis_func_grad(const Cell &c, const Node &n)
{
    ElementArray<Node> nodes = c.getNodes();
	double x[3];
	double y[3];
	// gradient of the basis function
	rMatrix grad(2,1);
	unsigned n_ind = 0;
	for(unsigned i = 0; i < 3; i++){
		if(n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}
	
	if(n_ind == 0){
		grad(0,0) = (y[1] - y[2]);
		grad(1,0) = - (x[1] - x[2]);
		grad /= ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2]));
	}
	else if(n_ind == 1){
		grad(0,0) = (y[0] - y[2]);
		grad(1,0) = - (x[0] - x[2]);
		grad /= ((x[1] - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y[1] - y[2]));
	}
	else if(n_ind == 2){
		grad(0,0) = (y[1] - y[0]);
		grad(1,0) = - (x[1] - x[0]);
		grad /= ((x[2] - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y[2] - y[0]));
	}
	else{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
	return grad;
}

void coords_from_barycentric(double *node_x, double *node_y, double *eta, double *x, double *y)
{
	*x = node_x[0] * eta[0] + node_x[1] * eta[1] + node_x[2] * eta[2];
	*y = node_y[0] * eta[0] + node_y[1] * eta[1] + node_y[2] * eta[2];
}

double integrate_over_triangle(const Cell &c, const Node n, double (*f)(double, double, const Cell&, const Node&))
{
	double res = 0.0;
	double w3 = 0.205950504760887;
	double w6 = 0.063691414286223;
	double eta3[3] = {0.124949503233232, 0.437525248383384, 0.437525248383384};
	double eta6[3] = {0.797112651860071, 0.165409927389841, 0.037477420750088};

	ElementArray<Node> nodes = c.getNodes();
	if(nodes.size() != 3){
		printf("Cell is not a triangle, has %lld nodes!\n", nodes.size());
		exit(1);
	}
	// Coordinates of triangle nodes
	double node_x[3], node_y[3];
	// Set them
	for(unsigned i = 0; i < 3; i++){
		double c[3];
		nodes[i].Centroid(c);
		node_x[i] = c[0];
		node_y[i] = c[1];
	}

	// Add contribution from all combinations in eta3
	double x, y, val;
	double eta[3];
	eta[0] = eta3[0];
	eta[1] = eta3[1];
	eta[2] = eta3[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	//printf("x = %e, y = %e, val = %e\n", x, y, val);
	res += w3 * val;
	eta[0] = eta3[1];
	eta[1] = eta3[2];
	eta[2] = eta3[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w3 * val*val;
	eta[0] = eta3[2];
	eta[1] = eta3[0];
	eta[2] = eta3[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w3 * val*val;


	// Add contribution from all combinations in eta6
	eta[0] = eta6[0];
	eta[1] = eta6[1];
	eta[2] = eta6[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[0];
	eta[1] = eta6[2];
	eta[2] = eta6[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[1];
	eta[1] = eta6[0];
	eta[2] = eta6[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[1];
	eta[1] = eta6[2];
	eta[2] = eta6[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[2];
	eta[1] = eta6[0];
	eta[2] = eta6[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[2];
	eta[1] = eta6[1];
	eta[2] = eta6[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;

	res *= c.Volume();
	return res;
}

double fphi(double x, double y, const Cell& c, const Node& n) {
	return source(x, y) * basis_func(c, n, x, y);
}

// [ [phi1,phi1], [phi1,phi2]...   ]
// [   ]
// [   ]
// energetic scalar product: [u,v] = (D*grad(u), grad(v)) = (grad(v))^T * D * grad(u)
void Problem::assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc)
{
     rMatrix D(2,2);
	 D(0,0) = c.RealArray(tagD)[0];
	 D(1,1) = c.RealArray(tagD)[1];
	 D(0,1) = D(1,0) = c.RealArray(tagD)[2];
	 
	 ElementArray<Node> nodes = c.getNodes();
	 
	 A_loc = rMatrix(3,3);
	 rhs_loc = rMatrix(3,1);
	 
	 for(unsigned i = 0; i < 3; i++){
		 rMatrix grad_phi_i = basis_func_grad(c, nodes[i]);
		 rhs_loc(i,0) = integrate_over_triangle(c, nodes[i], fphi);
		 for(unsigned j = 0; j < 3; j++){
			 rMatrix grad_phi_j = basis_func_grad(c, nodes[j]);
			double val = (grad_phi_j.Transpose() * (D * grad_phi_i))(0,0);
			A_loc(i,j) = val * c.Volume();
		 }
		 
		 // integrate RHS
	 }
	 
}

void Problem::assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs)
{
    /*
	// Cell loop
	// For each cell assemble local system
	// and incorporate it into global
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		rMatrix A_loc, rhs_loc;
		assembleLocalSystem(c, A_loc, rhs_loc);
		// Now A_loc is 3x3, rhs_loc is 3x1
		//  
		//
		//
		//

		ElementArray<Node> nodes = c.getNodes();
		unsigned glob_ind[3];
		for(unsigned loc_ind = 0; loc_ind < 3; loc_ind++){
			glob_ind[loc_ind] = static_cast<unsigned>(nodes[loc_ind].Integer(tagGlobInd));
		}
		
		for(unsigned loc_ind = 0; loc_ind < 3; loc_ind++){
			// Consider node with local index 'loc_ind'
			
			// Check if this is a Dirichlet node
			if(nodes[loc_ind].GetMarker(mrkDirNode))
				continue;
			
			for(unsigned j = 0; j < 3; j++){
				if(nodes[j].GetMarker(mrkDirNode)){
					rhs[glob_ind[loc_ind]] -= A_loc(loc_ind,j) * nodes[j].Real(tagBCval);
				}
				else
					A[glob_ind[loc_ind]][glob_ind[j]] += A_loc(loc_ind,j);
				
			}
			rhs[glob_ind[loc_ind]] += rhs_loc(loc_ind,0);
		}
	}
    */

    for(Mesh::iteratorFace face = m.BeginFace(); face != m.EndFace(); face++){
        Face f = face->getAsFace();

        //find n_f, x_f
        double n_f[3], x_f[3];
        f.UnitNormal(n_f);
        f.Barycenter(x_f);

        if(f.Boundary()){
            if (f.Integer(tagBCtype) == BC_DIR){
                Cell c = f.getCells()[0];
                double x_a[3];
                c.Barycenter(x_a);

                double diff_x, diff_y, diff_xy;
                diff_x = c.RealArray(tagD)[0];
                diff_y = c.RealArray(tagD)[1];
                diff_xy = c.RealArray(tagD)[2];

                double d_a[3];
                for (int i = 0; i < 3; i++){
                    d_a[i] = x_f[i] - x_a[i];
                }
                double kf;
                kf = (n_f[0] * (diff_x * d_a[0] + diff_xy * d_a[1]) + n_f[1] * (diff_xy * d_a[0] + diff_y * d_a[1]))
                        / (d_a[0] * d_a[0] + d_a[1] * d_a[1]);
                unsigned ind = static_cast<unsigned>(c.Integer(tagGlobInd));
                A[ind][ind] += kf * f.Area();
                rhs[ind] += kf * f.Area() * f.Real(tagBCval);

                f.Real(tagTrans) = kf;
            }
        }
        else{
            ElementArray<Cell> c = f.getCells();
            double x_a[3], x_b[3];
            c[0] = f.BackCell();
            c[1] = f.FrontCell();
            c[0].Barycenter(x_a);
            c[1].Barycenter(x_b);

            double diff_a_x, diff_a_y, diff_a_xy;
            double diff_b_x, diff_b_y, diff_b_xy;
            diff_a_x = c[0].RealArray(tagD)[0];
            diff_a_y = c[0].RealArray(tagD)[1];
            diff_a_xy = c[0].RealArray(tagD)[2];
            diff_b_x = c[1].RealArray(tagD)[0];
            diff_b_y = c[1].RealArray(tagD)[1];
            diff_b_xy = c[1].RealArray(tagD)[2];

            double d_a[3], d_b[3];
            for (int i = 0; i < 3; i++){
                d_a[i] = x_f[i] - x_a[i];
                d_b[i] = x_f[i] - x_b[i];
            }
            double kf_a, kf_b, tf;
            kf_a = (n_f[0] * (diff_a_x * d_a[0] + diff_a_xy * d_a[1]) + n_f[1] * (diff_a_xy * d_a[0] + diff_a_y * d_a[1]))
                   / (d_a[0] * d_a[0] + d_a[1] * d_a[1]);
            kf_b = (n_f[0] * (diff_b_x * d_b[0] + diff_b_xy * d_b[1]) + n_f[1] * (diff_b_xy * d_b[0] + diff_b_y * d_b[1]))
                   / (d_b[0] * d_b[0] + d_b[1] * d_b[1]);
            tf = (kf_a * kf_b) / (kf_a - kf_b);
            f.Real(tagTrans) = tf;

            unsigned ind0 = static_cast<unsigned>(c[0].Integer(tagGlobInd));
            unsigned ind1 = static_cast<unsigned>(c[1].Integer(tagGlobInd));
            A[ind0][ind0] += -tf * f.Area();
            A[ind0][ind1] += tf * f.Area();
            A[ind1][ind0] += tf * f.Area();
            A[ind1][ind1] += -tf * f.Area();
        }
    }



    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        unsigned ind = static_cast<unsigned>(c.Integer(tagGlobInd));
        rhs[ind] += c.Real(tagSource) * c.Volume();
    }
}

void Problem::run()
{
	// Matrix size
	//unsigned N = static_cast<unsigned>(m.NumberOfNodes()) - numDirNodes;
    unsigned N = static_cast<unsigned>(m.NumberOfCells());
	// Global matrix called 'stiffness matrix'
	Sparse::Matrix A;
	// Solution vector
	Sparse::Vector sol;
	// Right-hand side vector
	Sparse::Vector rhs;

	A.SetInterval(0, N);
	sol.SetInterval(0, N);
	rhs.SetInterval(0, N);

	assembleGlobalSystem(A, rhs);

	A.Save("A.mtx");
	rhs.Save("rhs.mtx");

	string solver_name = "inner_mptiluc";
	Solver S(solver_name);

	S.SetMatrix(A);
	bool solved = S.Solve(rhs, sol);
	if(!solved){
		printf("Linear solver failed: %s\n", S.GetReason().c_str());
		printf("Number of iterations: %d\n", S.Iterations());
		printf("Residual:             %e\n", S.Residual());
		exit(1);
	}

	for(Mesh::iteratorCell cell = m.BeginCell(); cell != m.EndCell(); cell++){
		Cell c = cell->getAsCell();
		if(c.GetMarker(mrkDirNode)){
			c.Real(tagConc) = c.Real(tagBCval);
			continue;
		}
		unsigned ind = static_cast<unsigned>(c.Integer(tagGlobInd));
		c.Real(tagConc) = sol[ind];
	}
	m.Save("res.vtk");
}

double Problem::C_norm() {
    double norm = 0.0;
    for(Mesh::iteratorCell cell = m.BeginCell(); cell != m.EndCell(); cell++){
        Cell c = cell->getAsCell();
        double diff = abs(c.Real(tagConc) - c.Real(tagConcAn));
        if (norm < diff){
            norm = diff;
        }
    }
    return norm;
}

double Problem::L2_norm() {
    double norm = 0.0;
    for(Mesh::iteratorCell cell = m.BeginCell(); cell != m.EndCell(); cell++){
        Cell c = cell->getAsCell();
        norm += abs(c.Real(tagConc) - c.Real(tagConcAn)) * c.Volume();
    }
    return norm;
}

int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.Load(argv[1]);

    main_mesh_diam(m);
	Problem P(m);
	P.initProblem();
	P.run();
    printf("C_norm = %lf\n", P.C_norm());
    printf("L2_norm = %lf\n", P.L2_norm());
	printf("Success\n");
	return 0;
}
