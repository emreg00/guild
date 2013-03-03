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

#ifndef EXCEPTIONS_HPP_29045684568920568920456590268
#define EXCEPTIONS_HPP_29045684568920568920456590268

#include <exception>
#include <stdexcept>

class TypeException: public std::exception
{
    virtual const char* what() const throw()
    {
	return "Type exception";
    }
} ; //typeEx;

class GenericError: public std::runtime_error
{
public:
    GenericError(const std::string& msg = ""):std::runtime_error(msg) {}
} ;

#endif // EXCEPTIONS_HPP_

