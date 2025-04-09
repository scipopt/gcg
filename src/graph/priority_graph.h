/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
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
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

struct Comparator{
    bool operator() (const pair<int, set<int> >& lhs, const pair<int, set<int> >& rhs){
        return lhs.second.size() > rhs.second.size()
                || (lhs.second.size() == rhs.second.size() && lhs.first > rhs.first);
    }
};

class priority_graph : public priority_queue<pair<int, set<int> >, vector<pair<int, set<int> > >, Comparator >
{
private:
    set<int> nodes;       // for optimization reasons
public:
    void addEdge(int node_i, int node_j){
        bool found1 = false;
        bool found2 = false;
        for(auto it = this->c.begin(); it < this->c.end(); ++it)
        {
            if(it->first == node_i){
                it->second.insert(node_j);
                found1 = true;
            }
            if(it->first == node_j){
                it->second.insert(node_i);
                found2 = true;
            }
            if(found1 && found2) break;
        }
        make_heap(this->c.begin(), this->c.end(), this->comp);
    }
    set<int> getNeighbors(int node){
        set<int> res;
        for(auto it = this->c.begin(); it < this->c.end(); ++it){
            if(it->first == node){
                return it->second;
            }
        }
        return res;
    }
    void addNode(int id){
        auto res = nodes.insert(id);
        if(res.second == true)
            this->push(pair<int, set<int> >(id, set<int>()));
    }
    bool removeNode(int node, vector<int>& removed){
        nodes.erase(node);
        bool res;
        auto it = this->c.begin();
        for(; it < this->c.end(); ++it)
        {
            if(it->first == node){
                break;
            }
        }
        if (it != this->c.end()) {
            this->c.erase(it);
            make_heap(this->c.begin(), this->c.end(), this->comp);
            removed.push_back(node);
            res = true;
        }
        else{
            cout << "failed to remove node " << node << endl;
            return false;
        }

        it = this->c.begin();
        for(; it < this->c.end(); ++it)
        {
            it->second.erase(node);
        }

        return res;
    }
    set<int> getNodes(){
       return nodes;
    }
};
