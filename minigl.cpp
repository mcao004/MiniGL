/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

// make structs
struct Vertex {
    vec4 pos;
    vec3 col;

    // constructors
    Vertex() : pos(vec<MGLfloat,4>(0,0,0,0)), col(vec<MGLfloat,3>(0,0,0)) {}
    Vertex(vec4 pos, vec3 col) : pos(pos), col(col) {}    
};

MGLfloat area(vec3 p1, vec3 p2, vec3 p3) {
    return cross(vec3(p2-p1),vec3(p3-p1)).magnitude() / 2;
}

struct Triangle {
    Vertex v1, v2, v3;

    // constructors
    Triangle() {}
    Triangle(Vertex v1, Vertex v2, Vertex v3) : v1(v1), v2(v2), v3(v3) {}
    
    // normal
    vec3 normal() {
        vec3 u1 = v2.pos.remove_last_dim()-v1.pos.remove_last_dim();
        vec3 u2 = v3.pos.remove_last_dim()-v1.pos.remove_last_dim();

        return vec<MGLfloat,3>(
                u1[1]*u2[2] - u1[2]*u2[1],
                u1[2]*u2[0] - u1[0]*u2[2],
                u1[0]*u2[1] - u1[1]*u2[0]);
    }

    // get Barycentric coordinates
    vec3 get_Bary(vec3 point) {
        MGLfloat a,b,c,area_tri; // alpha, beta, gamma and area of full triangle
        area_tri = area(v1.pos.remove_last_dim(),v2.pos.remove_last_dim(),v3.pos.remove_last_dim());
        a = area(point, v2.pos.remove_last_dim(),v3.pos.remove_last_dim()) / area_tri;
        b = area(v1.pos.remove_last_dim(), point,v3.pos.remove_last_dim()) / area_tri;
        c = 1.0f - a - b;

        return vec<MGLfloat,3>(a,b,c);
    }
};

MGLfloat x_min(Vertex v1, Vertex v2, Vertex v3) { return floor(min(min(v1.pos[0],v2.pos[0]),v3.pos[0])); }
MGLfloat x_max(Vertex v1, Vertex v2, Vertex v3) { return ceil(max(max(v1.pos[0],v2.pos[0]),v3.pos[0])); }
MGLfloat y_min(Vertex v1, Vertex v2, Vertex v3) { return floor(min(min(v1.pos[1],v2.pos[1]),v3.pos[1])); }
MGLfloat y_max(Vertex v1, Vertex v2, Vertex v3) { return ceil(max(max(v1.pos[1],v2.pos[1]),v3.pos[1])); }

// current color
vec3 curr_col;
MGLpoly_mode draw_mode;

// list of vertices and list of triangles
vector<Vertex> vert_list;
vector<Triangle> tri_list;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    Triangle curr_triangle;
    Vertex v1, v2, v3;
    MGLfloat xmin, xmax, ymin, ymax;

    cout << "width: " << width << endl;
    cout << "height: " << height << endl;

    for (unsigned tri = 0; tri < tri_list.size(); tri++) {
        curr_triangle = tri_list.at(tri);
        v1 = curr_triangle.v1;
        v2 = curr_triangle.v2;
        v3 = curr_triangle.v3;
        
        // transform the curr_triangle's vertices to screen space in v1,v2,v3
        v1.pos[0] = (width * (v1.pos[0] + 1)) / 2;
        v1.pos[1] = (height * (v1.pos[1] + 1)) / 2;
        v2.pos[0] = (width * (v2.pos[0] + 1)) / 2;
        v2.pos[1] = (height * (v2.pos[1] + 1)) / 2;
        v3.pos[0] = (width * (v3.pos[0] + 1)) / 2;
        v3.pos[1] = (height * (v3.pos[1] + 1)) / 2;

        // draw vertices to check
        cout << curr_triangle.v1.pos[0] << ", " << curr_triangle.v1.pos[1] << ", " << curr_triangle.v1.pos[2] << ", " << curr_triangle.v1.pos[3] << endl;
        cout << curr_triangle.v2.pos[0] << ", " << curr_triangle.v2.pos[1] << ", " << curr_triangle.v2.pos[2] << ", " << curr_triangle.v2.pos[3] << endl;
        cout << curr_triangle.v3.pos[0] << ", " << curr_triangle.v3.pos[1] << ", " << curr_triangle.v3.pos[2] << ", " << curr_triangle.v3.pos[3] << endl;

        // bounded box of triangle
        xmin = x_min(v1,v2,v3);
        xmax = x_max(v1,v2,v3);
        ymin = y_min(v1,v2,v3);
        ymax = y_max(v1,v2,v3);

        // vertice test
        data[(MGLint)(floor(v1.pos[0]) + width * floor(v1.pos[1]))] = Make_Pixel(255,255,255);
        data[(MGLint)(floor(v2.pos[0]) + width * floor(v2.pos[1]))] = Make_Pixel(255,255,255);
        data[(MGLint)(floor(v3.pos[0]) + width * floor(v3.pos[1]))] = Make_Pixel(255,255,255);
        
        // for each i and j in bounded box, move back to obj space and check if hits triangle
        MGLfloat denom;
        vec3 normal = curr_triangle.normal(); // normal of triangle
        vec3 proj_pt, bary, color; // projected point to its plane,  barycentric coords, and interpolated color
        for (MGLint i = xmin; i <= xmax; i++) {
            for (MGLint j = ymin; j <= ymax; j++) {
                data[width*j + i] = Make_Pixel(120,120,120);
                
                // i,j --> proj_pt in obj space
                proj_pt = vec<MGLfloat,3>( (((2 *(MGLfloat)i) + 1) / width) - 1,
                        (((2 *(MGLfloat)j) + 1) / height) - 1 , 0);

                // if inside current triangle
                denom = dot(normal,vec<MGLfloat,3>(0,0,3));
                if (denom < 1e-6) continue;

                proj_pt[2] = dot(vec3(v1.pos.remove_last_dim() - proj_pt), normal) / normal[2];

                //cout << "pix: " << i << ", " << j << " -----> " << 
                //    "proj_pt: " << proj_pt[0] << ", " << proj_pt[1] << ", " << proj_pt[2] << endl;
                
                // now we got a point on the same plane as the triangle
                // start barycentric stuff
                //proj_pt = curr_triangle.v1.pos.remove_last_dim();
                bary = curr_triangle.get_Bary(proj_pt);
                
                //cout << "Bary: " << bary[0] << ", " << bary[1] << ", " << bary[2] << endl;
                
                if (bary[0] >= 0 && bary[1] >= 0 && bary[2] >= 0) {
                    // cout << "int" << endl;
                    cout << "pix: " << i << ", " << j << " -----> " << 
                        "proj_pt: " << proj_pt[0] << ", " << proj_pt[1] << ", " << proj_pt[2] << endl;
                    cout << "Bary: " << bary[0] << ", " << bary[1] << ", " << bary[2] << endl;
                    
                    color = bary[0] * v1.col + bary[1] * v2.col + bary[2] * v3.col;
                    data[width*j + i] = Make_Pixel(color[0], color[1], color[2]);
                    data[width*j + i] = Make_Pixel(255,255,255);
                }
            }
        }


    }
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    draw_mode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
    unsigned i = 0;
    switch(draw_mode) {
        case MGL_TRIANGLES:
            while (i+2 < vert_list.size()) {
                tri_list.push_back(Triangle(vert_list.at(i), vert_list.at(i+1), vert_list.at(i+2)));
                i += 3;
            }
            break;
        case MGL_QUADS:
            while (i+3 < vert_list.size()) {
                tri_list.push_back(Triangle(vert_list.at(i), vert_list.at(i+1), vert_list.at(i+2)));
                tri_list.push_back(Triangle(vert_list.at(i), vert_list.at(i+2), vert_list.at(i+3))); 
                i += 4;
            }
            break;
    }
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    vec4 pos = vec<MGLfloat,4>(x,y,0,1);
    vec3 col = curr_col;
    vert_list.push_back(Vertex(pos, col));
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    vec4 pos = vec<MGLfloat,4>(x,y,z,1);
    vec3 col = curr_col; 
    vert_list.push_back(Vertex(pos,col));
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
    curr_col = vec<MGLfloat,3>(red,green,blue);
}
