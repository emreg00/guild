/*
    GUILD (Genes Underlying Inheritance Linked Disorders) implements several 
    graph based algorithms for scoring relevance of a node in the network in 
    terms of a phenotype using known associations in the node's neighborhood 
    for that phenotype. GUILD has been applied to the prioritization of genes 
    for several human disorders. 2011 - Emre Guney (Unviersitat Pompeu Fabra)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GLOBALS_H_
#define GLOBALS_H_
 
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

namespace fflow {

//typedef unsigned int uInt;

#define INT_INFINITY 9999999

typedef struct data_struct { 
	float reserve;
	float reserveUpdated;
	float capacity;
	data_struct() { reserve = 0.0; reserveUpdated = 0.0; capacity = 0.0; };
	data_struct(float valReserve, float valCapacity) { reserve = valReserve; reserveUpdated = valReserve; capacity = valCapacity; };
} Data;

}

#endif /*GLOBALS_H_*/

