#include <iostream>
#include <fstream>

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "Mesh.h"
#include "FVCA5Mesh.h"

using namespace std;

int main(int argc, char* argv[])
{ // nombre & tableaux de mots
  if (argc != 2) {
    cout << "Too much/few arguments... terminating." << endl;
    return EXIT_FAILURE;
  }

  // Reading mesh
  double scale = 100.;
  FVCA5Mesh Mh(argv[1] , scale);
  cout << "Opening " << argv[1] << " with scale=" << scale << std::endl;

  
  // 2D plot options
  sf::View view(sf::FloatRect(Mh.xmin, Mh.ymin, Mh.xmax, Mh.ymax));
  sf::RenderWindow window(sf::VideoMode(800,800), "Plot");
  window.setView(view); // focus the camera on the (0,100)x(0,100) region of the window
  window.setFramerateLimit(30); // FPS limiter
  sf::Color c_gray(200,200,200);
  sf::Color rtc[] = { sf::Color::Red, sf::Color::Blue };

  //
  sf::Event event;

  // redraw
  for (auto & M : Mh.Mailles) {
    M.shape.setFillColor(rtc[M.rt]); // oil sat in green , water sat in blue
    window.draw(M.shape); // draw in buffer
  }

  for (std::vector<Face>::iterator F=Mh.FacesStd.begin(); F != Mh.FacesStd.end(); F++) {
    sf::Vertex line[] = {
      sf::Vertex(sf::Vector2f(F->xa(0), Mh.ymax - F->xa(1)), c_gray),
      sf::Vertex(sf::Vector2f(F->xb(0), Mh.ymax - F->xb(1)), c_gray)
    };
    window.draw(line, 2, sf::Lines);
  }
  
  
  while (true) {

    // treat events 
    while (window.pollEvent(event)) { //

      if (event.type == sf::Event::MouseButtonReleased) {
	if (event.mouseButton.button == sf::Mouse::Left) {
	  // get the current mouse position in the window
	  sf::Vector2i pixelPos = sf::Mouse::getPosition(window);	 
	  sf::Vector2f worldPos = window.mapPixelToCoords(pixelPos);
	  
	  float px = worldPos.x;
	  float py = Mh.ymax - worldPos.y;
	  
	  Eigen::Vector2d xP;
	  xP << px, py;

	  // find corresponding cell
	  uint cell = 0;
	  double dist = 1e300;
	  for (auto & M : Mh.Mailles ) {	  
	    double s = (M.xT-xP).norm();	    
	    if ( s < dist ) {
	      cell = M.idx;
	      dist = s;
	    }
	  }
	  Mh.Mailles[cell].rt += 1;
	  Mh.Mailles[cell].rt %= 2;

	  
	}
      }
      
      if (event.type  == sf::Event::KeyPressed) {
	switch(event.key.code) {

	case sf::Keyboard::I :
	  {
	    for (auto &M : Mh.Mailles) {
	      M.rt += 1;
	      M.rt %= 2;
	    }
	  }
	  break;	  
	  
	case sf::Keyboard::S :
	  {
	    ifstream read(argv[1], ios::in);
	    ofstream write("out.tmp", ios::out); // read & write access

	    string marker;
	    read >> marker;
	    
	    while ( marker.compare("rocktype") != 0 and read.eof() != true) {
	      cout << "reading " << marker << "\n";
	      
	      write << marker << "\t";
	      cout << "writing " << marker << "\n\n";

	      // next
	      read >> marker;
	    }

	    read.close();
	    cout << "closing " << argv[1] << std::endl;

	    write << "rocktype" << "\t";
	    for (auto &M : Mh.Mailles) {
	      write << M.rt << "\t";
	    }
	    write.close();
	    cout << "closing " << "out.msh" << std::endl;

	    // rename 
	    
	  }
	  break;	  
	}
      }
    
    }
    
    // redraw
    for (auto & M : Mh.Mailles) {
      M.shape.setFillColor(rtc[M.rt]); // oil sat in green , water sat in blue
      window.draw(M.shape); // draw in buffer
    }

    for (std::vector<Face>::iterator F=Mh.FacesStd.begin(); F != Mh.FacesStd.end(); F++) {
      sf::Vertex line[] = {
	sf::Vertex(sf::Vector2f(F->xa(0), Mh.ymax - F->xa(1)), c_gray),
	sf::Vertex(sf::Vector2f(F->xb(0), Mh.ymax - F->xb(1)), c_gray)
      };
      window.draw(line, 2, sf::Lines);
    }
    window.display();
  }


  // terminating
  return EXIT_SUCCESS;
}
