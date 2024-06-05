/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Heur2opt.h
 * @brief  2-Optimum - combinatorial improvement heuristic for TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#pragma once
#include "util.hpp"
#include "instance.hpp"

#include "objscip/objscip.h"

/** C++ wrapper */
class MyScipHeur : public scip::ObjHeur
{

public:
    vector<vector<u_int32_t>> *var_index;
    vector<u_int32_t> *subgraph;
    OLCMInstance *instance;
    vector<u_int32_t> *index;
    vector<SCIP_Var *> *scipvars;
    vector<pair<u_int32_t, u_int32_t>>* var_vals;
    SCIP_Sol* sol;
    vector<u_int32_t> sortoffset;
    vector<u_int32_t> someorder;
    vector<vector<bool>> fixeddirection;

    /** default constructor */
    MyScipHeur(
        SCIP *scip,
        vector<vector<u_int32_t>> *var_index,
        vector<u_int32_t> *subgraph,
        OLCMInstance *instance,
        vector<u_int32_t> *index,
        vector<SCIP_Var *> *scipvars,
        vector<pair<u_int32_t, u_int32_t>>* var_vals);
    
    virtual ~MyScipHeur()
    {
    }
    
    virtual SCIP_DECL_HEURFREE(scip_free);
    
    virtual SCIP_DECL_HEURINIT(scip_init);
    
    virtual SCIP_DECL_HEUREXIT(scip_exit);
    
    virtual SCIP_DECL_HEURINITSOL(scip_initsol);
    
    virtual SCIP_DECL_HEUREXITSOL(scip_exitsol);
    
    virtual SCIP_DECL_HEUREXEC(scip_exec);
    
    // virtual SCIP_DECL_HEURCLONE(ObjCloneable *clone); /*lint !e665*/
    
    // virtual SCIP_DECL_HEURISCLONEABLE(iscloneable)
    // {
    //     return TRUE;
    // }
}; /*lint !e1712*/
