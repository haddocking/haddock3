#include <cstdio>
#include <cstring>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


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
  char seg;
};

vector<Residue> res;

int main(int argc, char *argv[]) {
  char buf[2000];

  if (argc < 3) {
    fprintf(stderr,"ERROR: Too few arguments\n");
    fprintf(stderr, "Usage: contact <cutoff> <pdb file>\n");
    return 1;
  }

  char *filename = argv[2];
  float cutoff = atof(argv[1]);

  if (cutoff < 0 || cutoff > 100) {
    fprintf(stderr,"ERROR: Cutoff out of range\n");
    fprintf(stderr, "Usage: contact <cutoff> <pdb file>\n");
    return 1;
  }

  FILE *fil = fopen(filename, "r");
  if (fil == NULL) {
    fprintf(stderr, "ERROR: PDB file %s does not exist\n", filename);
    return 1;
  }
  int currnr = -99999;
  char currseg = '#';
  set<int> nonconv;
  while (!feof(fil)) {
    char code[10];
    char atom[5];
    if (!fgets(buf, 2000, fil)) break;
    sscanf(buf, "%s %*d %s", code, atom);

    //LINE CHANGED BY K.Visscher
    //now test for both old and new style hydrogen
    //translation: only read if record == "ATOM" and if either the first letter in the atomname is NOT H or if not ( the first letter is a number and the second letter is H )
    if (!strncmp(code,"ATOM", 4) && ( atom[0] != 'H' && !(  isdigit(atom[0]) && atom[1] == 'H' )  )  ) {
      int nr = atoi(buf + 22);
      char seg = buf[72];
      if (nr != currnr || seg != currseg) {
        Residue r;
        r.nr = nr;
	r.seg = seg;
	res.push_back(r);
	currnr = r.nr;
	currseg = r.seg;
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

  if (!res.size()) {fprintf(stderr, "ERROR: PDB file %s contains no residues\n"); return 1;}

  double cutoffsq = cutoff * cutoff;

  for (int n = 0; n < res.size(); n++) {
    vector<Coor3f> &c1 = res[n].coor;
    char seg1 = res[n].seg;
    for (int nn = n + 1; nn < res.size(); nn++) {
      char seg2 = res[nn].seg;
      if (seg1 == seg2) continue;
      vector<Coor3f> &c2 = res[nn].coor;
      for (int i = 0; i < res[n].coor.size(); i++) {
        //printf ("%d %c %s\n", res[n].nr, res[n].seg, res[n].atom[i].c_str());
        for (int ii = 0; ii < res[nn].coor.size(); ii++) {
	  double currdissq =
	    (c1[i].x - c2[ii].x) * (c1[i].x - c2[ii].x) +
	    (c1[i].y - c2[ii].y) * (c1[i].y - c2[ii].y) +
	    (c1[i].z - c2[ii].z) * (c1[i].z - c2[ii].z);
	   if (currdissq < cutoffsq) {
	     printf ("%d %c %s %d %c %s %f\n", res[n].nr, res[n].seg, res[n].atom[i].c_str(), res[nn].nr, res[nn].seg, res[nn].atom[ii].c_str(), sqrt(currdissq));
	   }
        }
      }
    }
  }
}
