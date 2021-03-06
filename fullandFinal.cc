// Draw an Icosahedron
// ECE4893/8893 Project 4
// MAYANK PARASAR

#include <iostream>
#include <math.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>

using namespace std;

#define NFACE 20
#define NVERTEX 12

#define X .525731112119133606 
#define Z .850650808352039932

// Global variables are declared here
int WinH = 512;
int WinW = 512;
static GLfloat rotX = 0.0;
static GLfloat rotY = 0.0;
static GLfloat rotZ = 0.0;
bool rotate;


static GLuint maxDepth = 4;
int depth;


// These are the 12 vertices for the icosahedron
static GLfloat vdata[NVERTEX][3] = {    
   {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
   {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
   {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
};

// These are the 20 faces.  Each of the three entries for each 
// vertex gives the 3 vertices that make the face.
static GLint tindices[NFACE][3] = { 
   {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
   {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
   {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
   {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} 
    };

int testNumber; // Global variable indicating which test number is desired


// functions declarations:
void Test1(void);
void Test2(void);
void Test3(void);
void Test4(void);
void Test5(int);
void Test6(int);


// Normalize a vector of non-zero length
void normalize(GLfloat v[3])
{
  GLfloat d = sqrt (v[0] * v[0] +
                    v[1] * v[1] +
                    v[2] * v[2]);
  if (d == 0.0) return;
  v[0] /=d;
  v[1] /=d;
  v[2] /=d;
}
 
void drawTriangle(GLfloat* v1, GLfloat* v2, GLfloat* v3) 
{
  glBegin(GL_LINE_LOOP);
    glColor3f(1.0, 1.0, 1.0); //White lines
    glLineWidth(1.0);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);  
  glEnd();

  // Make color of face of triangle as function
  // of cordinates of vertex..

  GLfloat tmp_v1;
  GLfloat tmp_v2;
  GLfloat tmp_v3;
  tmp_v1 = (*v1 < 0) ? -(*v1) : *v1;
  tmp_v2 = (*v2 < 0) ? -(*v2) : *v2;
  tmp_v3 = (*v3 < 0) ? -(*v3) : *v3;

  glColor3f(tmp_v1, tmp_v2, tmp_v3);
  
  glBegin(GL_TRIANGLES);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
  glEnd();
  

} 

 
// recursively subdivide face 'depth' times
// and draw the resulting triangles
void subDivide(GLfloat* v1, GLfloat* v2, GLfloat* v3, int depth)
{
  if (depth == 0)
  {
    drawTriangle(v1, v2, v3); 
    return;
  }
  // Find the mid-point of each triangle side
  GLfloat v12[3];
  GLfloat v23[3];
  GLfloat v31[3];
  // calculate midpoints of each side 
  for (int i = 0; i < 3; ++i)
  {
    v12[i] = (v1[i] + v2[i])/2.0;
    v23[i] = (v2[i] + v3[i])/2.0;
    v31[i] = (v3[i] + v1[i])/2.0;
  }
  // excrude midpoints of to lie on unit sphere
  normalize(v12);
  normalize(v23);
  normalize(v31);

  // recursively subdivide new triangles
  subDivide(v1, v12, v31, depth-1);
  subDivide(v2, v23, v12, depth-1);
  subDivide(v3, v31, v23, depth-1);
  subDivide(v12, v23, v31, depth-1);

}

void my_display(void)
{
  static int pass;
  
  cout << "Displaying pass " << ++pass << endl;
  glEnable(GL_LINE_SMOOTH); // enable anti-aliasing
  // clear all
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // Clear the matrix
  glLoadIdentity();
  // Set the viewing transformation
  gluLookAt(0.0, 0.0, 1.0, // eye
            0.0, 0.0, 0.0, // at
            0.0, 1.0, 0.0); // up

  glTranslatef(250, 250, 250);
  glScalef(200, 200, 200);
  if (rotate) 
  {
    glRotatef(rotX, 1.0, 0.0, 0.0);
    glRotatef(rotY, 0.0, 1.0, 0.0);
    rotX += 1.0;
    rotY -= 1.0;  
  }

  // Draw icosahedron
  for (int i = 0; i < 20; i++) 
  {
    subDivide(vdata[tindices[i][0]],
              vdata[tindices[i][1]],
              vdata[tindices[i][2]], 
              depth/*Number of sub-division(depth)*/);
  }

  //Flush buffer
  glutSwapBuffers(); // If double buffering
}

void reshape(int w, int h)
{

  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, (GLdouble)w, (GLdouble)0.0, h, (GLdouble)w, (GLdouble)-w);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

}

void init()
{
  //select clearing (background) color
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_FLAT);
  // Enable depth buffer
  glEnable(GL_DEPTH_TEST);
}

void timer(int)
{
  glutPostRedisplay();
  glutTimerFunc(100.0, timer, 0);
}
// Test cases.  Fill in your code for each test case
void Test1()
{
  depth = 0;
  rotate = false;
}

void Test2()
{
  depth = 0;
  rotate = true;
}

void Test3()
{
  depth = 1;
  rotate = false;
}

void Test4()
{
  depth = 1;
  rotate = true;  
}

void Test5(int depth)
{
  ::depth = depth;
  rotate = false;
}

void Test6(int depth)
{
  ::depth = depth;
  rotate = true;
}


int main(int argc, char** argv)
{
  if (argc < 2)
    {
      std::cout << "Usage: icosahedron testnumber" << endl;
      exit(1);
    }
  // Set the global test number
  testNumber = atol(argv[1]);
  // Initialize glut  and create your window here
  // Set your glut callbacks here
  // Enter the glut main loop here
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(WinW, WinH);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("Icosahedron");

  // Initialize GL
  init();

  // set callbacks (aka Register functions here)
  glutDisplayFunc(my_display);
  glutReshapeFunc(reshape);
  glutTimerFunc(100.0, timer, 0);


  switch (testNumber) 
  {
    case 1: Test1();
            break;
    case 2: Test2();
            break;          
    case 3: Test3();
            break;
    case 4: Test4();
            break;
    case 5: depth = atoi(argv[2]);
            Test5(depth);
            break;          
    case 6: depth = atoi(argv[2]);
            Test6(depth);
            break;
    default: cout << "no maintching Test-number--- Exiting" << endl;
            exit(0); 
  }

  glutMainLoop(); 

  return 0;
}

