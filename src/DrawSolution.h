#include "Mesh.h"
#include "Data.h"
#include <SFML/Window.hpp>


sf::Color c_gray(220,220,220);
sf::Color c_black(0,0,0);

void DrawSolution(Mesh* Mh, Physics & P, Dat & X, sf::RenderWindow & window) {

   for (auto & M : Mh->Mailles) {
      R s = std::abs(X.S(M.idx));
      M.shape.setFillColor(sf::Color(120*s + (1-s)*255, (1-s)*255, (1-s)*255)); // 
      window.draw(M.shape); // draw in buffer
    }    

    
    // Boundary faces
    for (std::vector<BFace>::iterator F=Mh->BoundaryFaces.begin(); F != Mh->BoundaryFaces.end(); F++) {
      sf::Vertex line[] = {
	sf::Vertex(sf::Vector2f(F->xa(0), P.L - F->xa(1)), c_gray),
	sf::Vertex(sf::Vector2f(F->xb(0), P.b - F->xb(1)), c_gray)
      };
      window.draw(line, 2, sf::Lines);
    }	
    
    // Interior faces
    for (std::vector<Face>::iterator F=Mh->FacesStd.begin(); F != Mh->FacesStd.end(); F++) {
      sf::Vertex line[] = {
	sf::Vertex(sf::Vector2f(F->xa(0), P.L - F->xa(1)), c_gray),
	sf::Vertex(sf::Vector2f(F->xb(0), P.b - F->xb(1)), c_gray)
      };
      window.draw(line, 2, sf::Lines);
    }
   
    // Interfaces
    for (std::vector<IFace>::iterator F=Mh->Interfaces.begin(); F != Mh->Interfaces.end(); F++) {
      sf::Vertex line[] = {
	sf::Vertex(sf::Vector2f(F->xa(0), P.L - F->xa(1)), c_black), // interface en noir
	sf::Vertex(sf::Vector2f(F->xb(0), P.b - F->xb(1)), c_black)  // interface en noir
      };
      window.draw(line, 2, sf::Lines);
    }    
}
