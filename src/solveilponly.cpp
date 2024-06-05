
/*! \file tsp-subtour.cpp
  \brief Solves the TSP using dinamically generated subtour elimination constraints

  Solver for the traveling salesman problem using lazy constraints

*/
#include <OsiClpSolverInterface.hpp>
#include <CbcModel.hpp>
#include <OsiAuxInfo.hpp>
#include <CglCutGenerator.hpp>
#include <CoinTime.hpp>
#include "instance.hpp"
#include "preprocess.hpp"
#include "readdata.hpp"

// when separating fractional cuts
#define MIN_VIOLATION 1e-5

int verbose = 1;

int time_limit = 3600;

using namespace std;

typedef struct
{
  OLCMInstance *instance;
  int **x; // indexes of x variables
} TSPInfo;

// receives info of tsp instance and an integer solution, returns size of
// subtour and fills elements in els
int find_subtour(const TSPInfo *tspi, const double *s, int *els, int start);

// receives info of tsp instance and a fractional solution, returns the
// size of the subtour found and fills elements in els
int find_subtour_frac(const TSPInfo *tspi, const double *s, int *els, int start);

// for debugging
void print_sol(int n, const int **x, const double *s);

// check in a solution which is the output arc of node i
int out_arc(int n, const double *s, const int **x, int i);

class CglSubTour : public CglCutGenerator
{
public:
  CglSubTour(const TSPInfo *info);

  CglSubTour(const CglSubTour &rhs);

  virtual CglCutGenerator *clone() const;

  virtual void
  generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
               const CglTreeInfo info = CglTreeInfo());

  const TSPInfo *info;
};

int main()
{

  double start = CoinCpuTime();
  int **x; // indexes of x variables

  // reading distance matrix
  OLCMInstance instance = read_instance();
  reduce_twins(instance);
  compute_cuv(instance);
  int n = instance.working_nodes.size();
  int numpairs = (n * (n - 1)) / 2;

  x = new int *[n];
  x[0] = new int[n * n];
  for (int i = 1; (i < n); ++i)
    x[i] = x[i - 1] + n;

  int *idx = new int[numpairs];
  double *coef = new double[numpairs];
  double *lb = new double[numpairs];
  fill(lb, lb + numpairs, 0.0);
  double *ub = new double[numpairs];
  fill(ub, ub + numpairs, 1.0);
  double *obj = new double[numpairs];
  char name[256];
  vector<string> cnames;

  // creating x variables
  OsiClpSolverInterface *mip = new OsiClpSolverInterface();
  mip->messageHandler()->setLogLevel(0);
  int k = 0;

  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      sprintf(name, "x(%d,%d)", i, j);
      cnames.push_back(name);
      obj[k] = ((int)instance.CUV[i][j]) - ((int)instance.CUV[j][i]);

      instance.LB += instance.CUV[j][i];
      x[i][j] = k;
      ++k;
    }
  }
  CoinBigIndex *starts = new CoinBigIndex[numpairs + 1];
  fill(starts, starts + numpairs + 1, 0);
  mip->addCols(k, starts, NULL, NULL, lb, ub, obj);

  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      for (int k = j+1; k < n; k++) {
        idx[0] = x[i][j];
        idx[1] = x[j][k];
        idx[2] = x[i][k];
        coef[0] = 1.0;
        coef[1] = 1.0;
        coef[2] = -1.0;
        mip->addRow(3, idx, coef, 0, 1);
        break;
      }
    }
  }

  mip->setColNames(cnames, 0, cnames.size(), 0);
  for (int i = 0; (i < mip->getNumCols()); ++i)
    mip->setInteger(i);
  CbcModel model(*mip);

  // initial formulation is incomplete
  OsiBabSolver defaultC;
  defaultC.setSolverType(4);
  model.solver()->setAuxiliaryInfo(&defaultC);
  model.passInSolverCharacteristics(&defaultC);
  model.setKeepNamesPreproc(true);

  TSPInfo info;
  info.instance = &instance;
  info.x = x;

  CglSubTour cglst(&info);

  model.addCutGenerator(&cglst, 1, "subtour", true, true);

  model.setMaximumSeconds(time_limit);
  model.branchAndBound();

  cout << model.getBestPossibleObjValue() + instance.LB << endl;

  delete[] idx;
  delete[] coef;
  delete[] lb;
  delete[] ub;
  delete[] obj;

  delete mip;
  if (x)
  {
    delete[] x[0];
    delete[] x;
  }

  exit(0);
}

CglSubTour::CglSubTour(const TSPInfo *tspi)
{
  this->info = tspi;
}

CglSubTour::CglSubTour(const CglSubTour &rhs)
{
  this->info = rhs.info;
}

CglCutGenerator *CglSubTour::clone() const
{
  return new CglSubTour(this->info);
}

void CglSubTour::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
                              const CglTreeInfo tinfo)
{
  const TSPInfo *tspi = this->info;
  const double *s = si.getColSolution();
  const int **x = (const int **)tspi->x;
  int n = tspi->instance->working_nodes.size();
  int *idx = new int[3];
  double *coef = new double[3];

  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      for (int k = j + 1; k < n; k++)
      {
        idx[0] = x[i][j];
        idx[1] = x[j][k];
        idx[2] = x[i][k];
        coef[0] = 1.0;
        coef[1] = 1.0;
        coef[2] = -1.0;
        OsiRowCut cut;
        cut.setLb(-0.1);
        cut.setUb(1.1);
        cut.setRow(3, idx, coef);

        cs.insertIfNotDuplicate(cut);
      }
    }
  }

  idx[0] = 0;
  coef[0] = 1;

  delete[] idx;
  delete[] coef;
}

int find_subtour(const TSPInfo *tspi, const double *s, int *els, int start)
{
  return 0;
}

int find_subtour_frac(const TSPInfo *tspi, const double *s, int *els, int start)
{
  return 0;
}

void print_sol(int n, const int **x, const double *s)
{
  char col[256], str[256];
  printf("    ");
  for (int i = 0; (i < n); ++i)
  {
    sprintf(str, "%d", i);
    printf("%7s ", str);
  }
  printf("\n");
  for (int i = 0; (i < n); ++i)
  {
    sprintf(str, "%d", i);
    printf("%3s ", str);

    for (int j = 0; (j < n); ++j)
      printf("%.5f ", s[x[i][j]]);
    printf("\n");
  }
  printf("\n");
}

// check in a solution which is the output arc of node i
int out_arc(int n, const double *s, const int **x, int i)
{
  for (int j = 0; j < n; ++j)
  {
    if (s[x[i][j]] >= 0.99)
      return j;
  }
  return -1;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
