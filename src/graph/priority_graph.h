/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
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
