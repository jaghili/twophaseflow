#ifndef ELEMENT_H
#define ELEMENT_H

#include "defs.h"
#include <SFML/Graphics.hpp>

struct Element {
  uint idx; // indice de l'element
  
  iVec faces; // list of face global indices
  Vec xT; // center of the Element
  R vol; // volume of the element
  Eigen::MatrixXd sommets;
  sf::ConvexShape shape;
  
  uint rt;
};

#endif
