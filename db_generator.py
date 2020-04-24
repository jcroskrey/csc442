import networkx as nx
import sympy as sym
from functools import reduce
from copy import deepcopy
import pandas as pd
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, JSON, ForeignKey, create_engine, MetaData
from sqlalchemy.orm import sessionmaker
import csv
import json

Base = declarative_base()

engine = create_engine('sqlite:///poly.db', echo=False)

Session = sessionmaker(bind=engine, autoflush=False)
session = Session()

metadata = MetaData(engine)

sym.init_printing()
x = sym.Symbol('x')

class Ind(Base):
    __tablename__ = 'independence_poly'

    id = Column(Integer, primary_key=True)
    edge_set = Column(JSON, unique=True)
    polynomial = Column(String)
    coefficients = Column(String)
    degree = Column(Integer)

    def __repr__(self):
        return "<Ind(edge_set='%s', polynomial='%s', coefficients='%s', \
                degree='%s')>" % (
                    self.edge_set, self.polynomial, self.coefficients, self.degree
                )


Base.metadata.create_all(engine)


def delete_node(T):
    """Deletes largest numerical node in T.
    
    Returns T.
    """
    size = len(T)
    T.remove_node(size-1)
    return T      

def delete_neighborhood(T):
    """Deletes neighborhood of largest numerical node in T.
    
    Returns T.
    """
    size = len(T.edges)
    neighbors = [n for n in T[size]]
    
    for neighbor in neighbors:
        T.remove_node(neighbor)
    
    T.remove_node(size)
    
    return T

def components(comp):
    """Takes connected components of the deleted neighborhood tree and passes them to 
    
    minus_neighborhood_polys.
    """
    S = [comp.subgraph(c) for c in nx.connected_components(comp)]
    return S

def polynomial_of_tree(T):
    """ Matches tree to independence polynomial within our dictionary.
    
    If the tree consists of one or two nodes, then it automatically assigns it's polynomial.
    """
    if len(T.nodes) == 1:
        f = x + 1
        return f

    if len(T.nodes) == 2:
        f = 2*x + 1
        return f

    if len(T.nodes) >= 3:
        length = len(T.nodes)

        for tree in nx.nonisomorphic_trees(length):
            if nx.is_isomorphic(tree, T):
                T = tree

        for row in session.query(Ind).filter_by(degree=length):
            if row.edge_set == json.dumps(list(T.edges)):
                f = sym.sympify(row.polynomial)
                return f


def multiply_polys(poly_list):
    """Multiplies list of polynomials together.
    """
    final_poly = reduce(lambda x, y: x*y, poly_list)
    
    return final_poly.expand()

def grab_coeffs(poly):
    p = sym.Poly(poly, x)
    c = p.coeffs()
    return c

def populate_independence_db(index):
    """The master plan that calls the functions and keeps track of different variables etc.
    """
    global db_index

    for tree in nx.nonisomorphic_trees(index):
        # First, create a copy of the object to work with.
        tree1 = tree 
        tree2 = deepcopy(tree)
        original_tree = deepcopy(tree2) # Keep this copy for when we add the tree to the dictionary.
        
        # Then deal with the first tree (tree1) where we delete the node. 
        delete_node(tree1)
        final_poly_tree1 = polynomial_of_tree(tree1)
        coefficients_poly1 = grab_coeffs(final_poly_tree1)
        # Note that the coefficients are reversed. They start at the highest degree of x.
        # This is because this is python and sympy's default way of handling polynomials.
        coefficients_poly1.insert(0, 0)
        
        # Lastly, delete the neighborhood of the node from tree2 and find polynomial.
        
        delete_neighborhood(tree2)
        components_tree2 = components(tree2)
        
        polys_of_tree2 = []
        for subgraph in components_tree2:
            polys_of_tree2.append(polynomial_of_tree(subgraph))
        
        # Multiply polynomials together    
        final_poly_tree2 = reduce(lambda x, y: x*y, polys_of_tree2)
        
        # Grab coefficients and add zero to the end of the list.
        coefficients_poly2 = grab_coeffs(final_poly_tree2)
        # Note that the coefficients are reversed. They start at the highest degree of x.
        # This is because this is python's and sympy's default way of handling polynomials.
        coefficients_poly2.append(0)
        
        while len(coefficients_poly1) > len(coefficients_poly2):
            coefficients_poly2.insert(0, 0)
        
        final_coeffs = [sum(i) for i in zip(coefficients_poly1, coefficients_poly2)]
        final_poly = sum(co*x**i for i, co in enumerate(reversed(final_coeffs)))
        
        # Add the data to our Ind database.
        new_entry = Ind(edge_set=json.dumps(list(original_tree.edges)), 
                        polynomial=str(final_poly), # To convert to sympy.core.add.Add use sym.sympify(poly)
                        coefficients=str(final_coeffs), # To convert back to list use eval(coeffs)
                        degree=len(original_tree)) 

        session.add(new_entry)

    return None

def main():
    """ Populate the dictionary.
    """
    for degree in range(8, 13):
        populate_independence_db(degree)
        session.commit()

    return None


if __name__ == "__main__":
    main()