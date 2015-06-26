#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

// use bit pattern to enumerate all possibilities

typedef struct{
  char Type[2];
  int Z;
  double R[3];
}Atom;

typedef struct{
  char Name[20];
  char Cell[50];
  int Scale[3], N, n, Ntype, *Zs, *Ns;
  Atom *Atoms;
}Crystal;

int type2Z(char*);
char* Z2type(int);
int Z2v(int);
void loadCrystal(char*, Crystal*);
void sortCrystal(Crystal*);
void shuffleCrystal(int, Crystal*);
void printCrystal(Crystal);
void printCPMD(Crystal);
void printCPMDFile(Crystal);

int main(int argn, char* argv[]){
  int i;
  Crystal cylInp;

  if(argn != 2){
    printf("expecting ONE argument for total number of atoms\n");
    return 1;
  }

  loadCrystal(argv[1], &cylInp);
  sortCrystal(&cylInp);
  //printCrystal(cylInp);
  //printCPMD(cylInp);
  for(i=0;i<pow(2,cylInp.n);i++){
    shuffleCrystal(i,&cylInp);
    sortCrystal(&cylInp);
    //printCrystal(cylInp);
    printCPMDFile(cylInp);
  }
  //printCrystal(cylInp);
}

int type2Z(char *type){
  if(strcmp(type,"Al")==0){
    return 13;
  }else if(strcmp(type, "Ga")==0){
    return 31;
  }else if(strcmp(type, "As")==0){
    return 33;
  }else{
    printf("ERROR: %s is not yet coded\n",type);
    return 0;
  }
}

char* Z2type(int Z){
  if(Z==13){
    return "Al";
  }else if(Z==31){
    return "Ga";
  }else if(Z==33){
    return "As";
  }else{
    printf("ERROR: Z=%d is not yet coded\n",Z);
    return 0;
  }
}

int Z2v(int Z){
  if(Z==13){
    return 3;
  }else if(Z==31){
    return 3;
  }else if(Z==33){
    return 5;
  }else{
    printf("ERROR: Z=%d is not yet coded\n",Z);
    return 0;
  }
}

void loadCrystal(char inp[20], Crystal *cylPtr){
  FILE *fr;
  size_t len=0;
  ssize_t read;
  char *string = (char *) malloc(80 * sizeof(char));
  int *tmp;
  int N, i, j, Z;

  fr = fopen(inp, "r");
  fscanf(fr,"%d",&N);
  read = getline(&string, &len, fr);
  read = getline(&string, &len, fr);
  strcpy(cylPtr->Cell,string);
  for(i=0;i<3;i++){
    fscanf(fr, "%d", &cylPtr->Scale[i]);
  }
  read = getline(&string, &len, fr);
  
  cylPtr->Atoms = malloc(N * sizeof(Atom));
  cylPtr->N = N;
  cylPtr->n = N/2;
  cylPtr->Ntype = 0;
  strcpy(cylPtr->Name, inp);

  for(i=0;i<N;i++){
    fscanf(fr, "%s %lf %lf %lf", 
      string, 
      &cylPtr->Atoms[i].R[0], 
      &cylPtr->Atoms[i].R[1], 
      &cylPtr->Atoms[i].R[2]);
    Z = type2Z(string);
    strcpy(cylPtr->Atoms[i].Type, string);
    cylPtr->Atoms[i].Z = Z;
  }
}

void sortCrystal(Crystal *sortCyl){
  int i, j, Ntype=1, *Zs, *Ns;
  Atom swap;

  void sort(Crystal *sortCyl, int k){
    int i;
    for(i=1;i<sortCyl->N;i++){
      for(j=i;j<sortCyl->N;j++){
        if(sortCyl->Atoms[j].Z == sortCyl->Atoms[i].Z){
          if(sortCyl->Atoms[j].R[k] < sortCyl->Atoms[i].R[k]){
            swap = sortCyl->Atoms[j];
            sortCyl->Atoms[j] = sortCyl->Atoms[i];
            sortCyl->Atoms[i] = swap;
          }
        }
      }
    }
  }

  for(i=0;i<sortCyl->N;i++){
    for(j=i;j<sortCyl->N;j++){
      if(sortCyl->Atoms[j].Z < sortCyl->Atoms[i].Z){
        swap = sortCyl->Atoms[j];
        sortCyl->Atoms[j] = sortCyl->Atoms[i];
        sortCyl->Atoms[i] = swap;
      }
    }
  }
  for(i=2;i>=0;i--) sort(sortCyl,i);

  for(i=1;i<sortCyl->N;i++){
    if(sortCyl->Atoms[i].Z != sortCyl->Atoms[i-1].Z){
      Ntype++;
    }
  }
  Zs = (int*) malloc(Ntype * sizeof(int));
  Ns = (int*) malloc(Ntype * sizeof(int));
  Zs[0] = sortCyl->Atoms[0].Z;
  for(i=0;i<Ntype;i++) Ns[i]=1;
  j=0;
  for(i=1;i<sortCyl->N;i++){
    if(sortCyl->Atoms[i].Z != sortCyl->Atoms[i-1].Z){
      Zs[++j] = sortCyl->Atoms[i].Z;
    }else{
      Ns[j]++;
    }
  }
  sortCyl->Zs = Zs;
  sortCyl->Ns = Ns;
  sortCyl->Ntype = Ntype;
}

void shuffleCrystal(int index, Crystal *cylShuffle){
  int i, N;

  for(i=cylShuffle->n-1; i>=0; i--){
    if(index&(1<<i)){
      // 1 denotes Al
      putchar('1');
      strcpy(cylShuffle->Atoms[i].Type, "Al");
      cylShuffle->Atoms[i].Z = 13;
    }else{// if(i<X){
      // 0 denotes Ga
      putchar('0');
      strcpy(cylShuffle->Atoms[i].Type, "Ga");
      cylShuffle->Atoms[i].Z = 31;
    }
  }
  sprintf(cylShuffle->Name, "v%04d%s",index+1, ".inp");
  printf("\n");
}

void printCrystal(Crystal cyl2Print){
  int i;
  printf("crystal file name: %s\n", cyl2Print.Name);
  printf("number of atoms: %d\n", cyl2Print.N);
  printf("number of atom types: %d\n", cyl2Print.Ntype);
  printf("atom list: ");
  for(i=0;i<cyl2Print.Ntype;i++){
    printf("%s ",Z2type(cyl2Print.Zs[i]));
  }
  printf("\n");
  for(i=0;i<cyl2Print.N;i++){
    printf("%2s %d % 5.3lf % 5.3lf % 5.3lf\n",
      cyl2Print.Atoms[i].Type, 
      cyl2Print.Atoms[i].Z,
      cyl2Print.Atoms[i].R[0], 
      cyl2Print.Atoms[i].R[1], 
      cyl2Print.Atoms[i].R[2]);
  }
}

void printCPMD(Crystal cyl2Print){
  int i, j, J0=0;

  //FILE *fr;

  //fr = fopen(cyl2Print.Name, "w");

  printf("&CPMD\n");
  printf(" OPTIMIXE WAVEFUNCTION\n");
  printf(" MIRROR\n");
  printf("&END\n");
  printf("\n");
  printf("&SYSTEM\n");
  printf(" SYMMETRY\n");
  printf("  ORTHORHMBIC\n");
  printf(" CELL ABSOLUTE\n");
  printf("  %s",cyl2Print.Cell);
  printf(" SCALE SX=%d SY=%d SZ=%d\n",
         cyl2Print.Scale[0],
         cyl2Print.Scale[1],
         cyl2Print.Scale[2]);
  printf(" ANGSTROM\n");
  printf(" CUTOFF\n");
  printf("  100\n");
  printf(" KPOINTS MONKHORST-PACK\n");
  printf("  1 1 1\n");
  printf("&END\n");
  printf("\n");
  printf("&ATOMS\n");
  for(i=0;i<cyl2Print.Ntype;i++){
    printf("*%s_q%d_pbe.psp FRAC\n",
           Z2type(cyl2Print.Zs[i]),
           Z2v(cyl2Print.Zs[i]));
    printf(" LMAX=F\n  %d\n", cyl2Print.Ns[i]);
    for(j=0;j<cyl2Print.Ns[i];j++){
      printf("   % 5.3lf % 5.3lf % 5.3lf\n",
        cyl2Print.Atoms[j+J0].R[0], 
        cyl2Print.Atoms[j+J0].R[1], 
        cyl2Print.Atoms[j+J0].R[2]);
    }
    printf("\n");
    J0 += cyl2Print.Ns[i];
  }
  printf("&END\n");
}

void printCPMDFile(Crystal cyl2Print){
  int i, j, J0=0;

  FILE *fr;

  fr = fopen(cyl2Print.Name, "w");

  fprintf(fr, "&CPMD\n");
  fprintf(fr, " OPTIMIXE WAVEFUNCTION\n");
  fprintf(fr, " MIRROR\n");
  fprintf(fr, "&END\n");
  fprintf(fr, "\n");
  fprintf(fr, "&DFT\n");
  fprintf(fr, " FUNCTIONAL PBE\n");
  fprintf(fr, "&END\n");
  fprintf(fr, "\n");
  fprintf(fr, "&SYSTEM\n");
  fprintf(fr, " SYMMETRY\n");
  fprintf(fr, "  ORTHORHMBIC\n");
  fprintf(fr, " CELL ABSOLUTE\n");
  fprintf(fr, "  %s",cyl2Print.Cell);
  fprintf(fr, " SCALE SX=%d SY=%d SZ=%d\n",
         cyl2Print.Scale[0],
         cyl2Print.Scale[1],
         cyl2Print.Scale[2]);
  fprintf(fr, " ANGSTROM\n");
  fprintf(fr, " CUTOFF\n");
  fprintf(fr, "  100\n");
  fprintf(fr, " KPOINTS MONKHORST-PACK\n");
  fprintf(fr, "  1 1 1\n");
  fprintf(fr, "&END\n");
  fprintf(fr, "\n");
  fprintf(fr, "&ATOMS\n");
  for(i=0;i<cyl2Print.Ntype;i++){
    fprintf(fr, "*%s_q%d_pbe.psp FRAC\n",
           Z2type(cyl2Print.Zs[i]),
           Z2v(cyl2Print.Zs[i]));
    fprintf(fr, " LMAX=F\n  %d\n", cyl2Print.Ns[i]);
    for(j=0;j<cyl2Print.Ns[i];j++){
      fprintf(fr, "   % 5.3lf % 5.3lf % 5.3lf\n",
        cyl2Print.Atoms[j+J0].R[0], 
        cyl2Print.Atoms[j+J0].R[1], 
        cyl2Print.Atoms[j+J0].R[2]);
    }
    fprintf(fr, "\n");
    J0 += cyl2Print.Ns[i];
  }
  fprintf(fr, "&END\n");
  fclose(fr);
}
