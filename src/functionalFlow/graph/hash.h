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

#ifndef HASH_H_
#define HASH_H_

#include "globals.h"
#include <unordered_map> //#include "hash_map.h"
#include <cstring>

namespace fflow {

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

struct isEqualString
{
  bool operator()(string s1, string s2) const
  {
    return s1.compare(s2) == 0;
  }
};

/*
class stringhasher : public stdext::hash_compare <std::string>
{
public:
    // Inspired by the java.lang.String.hashCode() algorithm somewhat processor cache-friendly), @param The string to be hashed, @return The hash value of s
    size_t operator() (const std::string& s) const {
        size_t h = 0;
        std::string::const_iterator p, p_end;
        for(p = s.begin(), p_end = s.end(); p != p_end; ++p) {
            h = 31 * h + (*p);
        }
        return h;
    }
    // @param s1 The first string, @param s2 The second string @return true if the first string comes before the second in lexicographical order
    bool operator() (const std::string& s1, const std::string& s2) const {
        return s1 < s2;
    }
};

typedef stdext::hash_map<std::string, std::string, stringhasher> HASH_S_S;

    hm[std::string("novembro")] = std::string("November");
    it = hm.find(std::string("março"));
    if(it != hm.end())
        std::cout << "The value corresponding to the key 'março' is " << it->second << std::endl;
    else
        std::cout << "The value corresponding to the key 'março' was not found" << std::endl;
*/

}

#endif /*HASH_H_*/

