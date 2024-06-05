#ifndef __SCIP_CONS_StrongTriangles_H__
#define __SCIP_CONS_StrongTriangles_H__

#include "instance.hpp"
#include "scip/scip.h"
#include "objscip/objscip.h"

class ConshdlrStrongTriangles : public scip::ObjConshdlr
{
public:
   vector<vector<bool>> isintersect;
   vector<vector<u_int32_t>> intersecting;
   vector<vector<u_int32_t>> *var_index;
   vector<u_int32_t> *subgraph;
   OLCMInstance *instance;
   vector<u_int32_t> *index;
   vector<SCIP_Var *> *scipvars;
   SCIP_Cons* dummycons;
   int stage;
   ConshdlrStrongTriangles(SCIP *scip,
                           vector<vector<u_int32_t>> *var_index,
                           vector<u_int32_t> *subgraph,
                           OLCMInstance *instance,
                           vector<u_int32_t> *index,
                           vector<SCIP_Var *> *scipvars) : ObjConshdlr(scip, "Strong triangles constraint handler", "Strong triangles constraint handler", 100, -100, -100, 1, 1, 100, -1, false, true, false, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_MEDIUM)
   {
      this->subgraph = subgraph;
      this->instance = instance;
      this->index = index;
      this->scipvars = scipvars;
      this->var_index = var_index;
      stage = 0;
      isintersect = vector<vector<bool>>(subgraph->size(), vector<bool>(subgraph->size(), false));
      intersecting = vector<vector<u_int32_t>>(subgraph->size());
      for (u_int32_t i = 0; i < subgraph->size(); i++)
      {
         auto u = (*subgraph)[i];

         for (u_int32_t j = 0; j < subgraph->size(); j++)
         {
            auto v = (*subgraph)[j];
            auto sa = instance->al[u][0];
            auto ea = instance->al[u][instance->al[u].size() - 1];
            auto sb = instance->al[v][0];
            auto eb = instance->al[v][instance->al[v].size() - 1];
            if (sa == eb || sb == ea)
               continue;
            // isintersect[i][j] = true;
            if ((sa <= sb && sb <= ea) || (sb <= sa && sa <= eb))
            {
               isintersect[i][j] = true;
            }
         }
      }
      for (u_int32_t i = 0; i < subgraph->size(); i++)
      {
         auto u = (*subgraph)[i];

         for (u_int32_t j = 0; j < subgraph->size(); j++)
         {
            if (i != j)
            {
               if (isintersect[i][j])
                  intersecting[i].push_back(j);
            }
         }
      }
   }

   virtual ~ConshdlrStrongTriangles()
   {
   }

   virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

   virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);

   virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

   // virtual SCIP_DECL_CONSENFORELAX(scip_enforelax);

   virtual SCIP_DECL_CONSENFOPS(scip_enfops);

   virtual SCIP_DECL_CONSCHECK(scip_check);

   virtual SCIP_DECL_CONSLOCK(scip_lock);

   virtual SCIP_DECL_CONSTRANS(scip_trans) {
      assert(false);
      return SCIP_OKAY;
   }

   // virtual SCIP_DECL_CONSPROP(scip_prop);

   // virtual SCIP_DECL_CONSINITLP(init);

   // void propagatetriangle(SCIP* scip, u_int32_t i, u_int32_t j, u_int32_t k, SCIP_RESULT* result, SCIP_Bool* tightened, SCIP_Bool* infeasible, SCIP_CONSHDLR* hdlr);

   bool is_constr_consistent_integer(SCIP* scip, SCIP_SOL* sol, const int i, const int j, const int k);

   bool is_consistent_integer(SCIP *scip, SCIP_SOL *sol);

   bool separate_if_efficatious(SCIP *scip, SCIP_SOL *sol, const int i, const int j, const int k, unsigned int enforce, SCIP_CONSHDLR* hdlr, double minviolation, SCIP_RESULT* result);

   bool storecut(SCIP *scip, SCIP_SOL *sol, const int i, const int j, const int k, unsigned int enforce, SCIP_CONSHDLR* hdlr, double minviolation, map<double, vector<SCIP_ROW *>>& cuts);

   SCIP_RETCODE sepaStrongTriangles(SCIP* scip, SCIP_SOL* sol, SCIP_Bool enforce, SCIP_RESULT* result, SCIP_Conshdlr* hdlr);

   

   /**@addtogroup CONSHDLRS
    *
    * @{
    *
    * @name StrongTriangles Constraints
    *
    * @{
    */
};

/** creates and captures a StrongTriangles constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE
SCIPcreateConsBasicStrongTriangles(
    SCIP *scip,       /**< SCIP data structure */
    SCIP_CONS **cons, /**< pointer to hold the created constraint */
    const char *name, /**< name of constraint */
    int nvars,        /**< number of variables in the constraint */
    SCIP_VAR **vars,  /**< array with variables of constraint entries */
    SCIP_Real *coefs, /**< array with coefficients of constraint entries */
    SCIP_Real lhs,    /**< left hand side of constraint */
    SCIP_Real rhs     /**< right hand side of constraint */
);

#endif
