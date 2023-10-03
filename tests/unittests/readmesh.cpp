#include <fstream>

//-----------------------------------------------------------------------------
void read_off(const char* fname, int &nv, int &nt, double* &ver, int* &tri)
{
    char buffer [129];
    int dummy;
    std::ifstream file(fname);
    file >> buffer;
    file >> nv >> nt >> dummy;
    ver = new double [3*nv];
    for (int i=0;i<3*nv;i++) file >> ver[i];
    tri = new int [3*nt];
    for (int i=0;i<nt;i++) {
        file >> dummy >> tri[3*i] >> tri[3*i+1] >> tri[3*i+2];
    }
    file.close();
}

//-----------------------------------------------------------------------------
void read_orient_txt(const char* fname, int nt, double* &orient)
{
    orient = new double [3*nt];
    std::ifstream file(fname);
    for (int i=0;i<3*nt;i++) file >> orient[i];
    file.close();
}

//-----------------------------------------------------------------------------
void dummy_orient(int nt, double* &orient)
{
    orient = new double [3*nt];
    for (int i=0;i<nt;i++) {
        orient[3*i] = 1;
        orient[3*i+1] = orient[3*i+2] = 0;
    }
}