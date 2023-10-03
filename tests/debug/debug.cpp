void write_streamline_collection(const char* fname, std::vector<StreamlineContainer*> collection)
{
    FILE* fid = fopen(fname, "w");
    for (unsigned int k=0;k<collection.size();k++) {
        for (unsigned int i=0;i<collection[k]->n;i++) {
            fprintf(fid, "%.5f %.5f %.5f  ", collection[k]->points[3*i], 
                                             collection[k]->points[3*i+1], 
                                             collection[k]->points[3*i+2]);
        }
        fprintf(fid, "\n");
    }
    fclose(fid);
}

void randomize_orient(int nt, double* orient, bool sign, bool all)
{
    if (all)
        for (int i=0;i<nt;i++) {
            orient[3*i]   = rand()/((double)RAND_MAX);
            orient[3*i+1] = rand()/((double)RAND_MAX);
            orient[3*i+2] = 0;
        }
    if (sign)
        for (int i=0;i<nt;i++) {
            if (rand() < RAND_MAX / 2) {
                orient[3*i] = -orient[3*i];
                orient[3*i+1] = -orient[3*i+1];
                orient[3*i+2] = -orient[3*i+2];
            }
        }
}
