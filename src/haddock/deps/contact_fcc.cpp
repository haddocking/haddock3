#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

struct Coor3f {
  float x;
  float y;
  float z;
};

struct Residue {
  int nr;
  vector<Coor3f> coor;
  vector<string> atom;
  int seg;
};

vector<Residue> res;

bool seg_sorter (Residue res_a, Residue res_b) {
  int segA = res_a.seg;
  int segB = res_b.seg;
  return (segA < segB);
};

int main(int argc, char *argv[]) {
  char buf[2000];

  if (argc < 3) {
    fprintf(stderr,"ERROR: Too few arguments\n");
    fprintf(stderr, "Usage: contact <pdb file> <cutoff>\n");
    return 1;
  }

  char *filename = argv[1];
  float cutoff = atof(argv[2]);

  if (cutoff < 0 || cutoff > 100) {
    fprintf(stderr,"ERROR: Cutoff out of range\n");
    fprintf(stderr, "Usage: contact <pdb file> <cutoff>\n");
    return 1;
  }

  FILE *fil = fopen(filename, "r");
  if (fil == NULL) {
    fprintf(stderr, "ERROR: PDB file %s does not exist\n", filename);
    return 1;
  }
  
  int currnr = -99999;
  char currseg;
  int segid = 0;

  set<int> nonconv;
  while (!feof(fil)) {
    char code[10];
    char atom[5];
    if (!fgets(buf, 2000, fil)) break;
    sscanf(buf, "%s %*d %s", code, atom);

    // Ignore HETATM and hydrogens
    if (!strncmp(code,"ATOM", 4) && ( atom[0] != 'H' && !(  isdigit(atom[0]) && atom[1] == 'H' )  )  ) {
      int nr = atoi(buf + 22);
      char seg = buf[72];
      if (seg != currseg) {
        currseg = seg;
        segid++;
      }
      if (nr != currnr) {
          Residue r;
          r.nr = nr+10000;
    	    r.seg = segid;
  	      res.push_back(r);
  	      currnr = r.nr;
      }
      Residue &rcurr = res[res.size() -1];
      Coor3f ccoor;
      ccoor.x = atof(buf+27);
      ccoor.y = atof(buf+38);
      ccoor.z = atof(buf+46);
      rcurr.coor.push_back(ccoor);
      string atom2(atom);
      rcurr.atom.push_back(atom2);
    }
  }

  if (!res.size()) {fprintf(stderr, "ERROR: PDB file %s contains no residues\n", filename); return 1;}

  // Sort the residues by segment to avoid random chain ordering problems
  sort (res.begin(), res.end(), seg_sorter);

  double cutoffsq = cutoff * cutoff;

  for (int n = 0; n < res.size(); n++) {
    vector<Coor3f> &c1 = res[n].coor;
    int seg1 = res[n].seg;
    for (int nn = n + 1; nn < res.size(); nn++) {
      int seg2 = res[nn].seg;
      if (seg1 == seg2) continue;
      vector<Coor3f> &c2 = res[nn].coor;
      for (int i = 0; i < res[n].coor.size(); i++) {
        for (int ii = 0; ii < res[nn].coor.size(); ii++) {
	  double currdissq =
	    (c1[i].x - c2[ii].x) * (c1[i].x - c2[ii].x) +
	    (c1[i].y - c2[ii].y) * (c1[i].y - c2[ii].y) +
	    (c1[i].z - c2[ii].z) * (c1[i].z - c2[ii].z);
	   if (currdissq < cutoffsq) {
	     printf ("%d%d%d%d\n", res[n].nr, res[n].seg, res[nn].nr, res[nn].seg);
	   }
        }
      }
    }
  }
  fclose(fil);
}
