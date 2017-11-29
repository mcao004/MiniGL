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

#include <stack>

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
        vec4 v1_pos = this->v1.pos / this->v1.pos[3];
        vec4 v2_pos = this->v2.pos / this->v2.pos[3];
        vec4 v3_pos = this->v3.pos / this->v3.pos[3];
        
        vec3 u1 = v2_pos.remove_last_dim()-v1_pos.remove_last_dim();
        vec3 u2 = v3_pos.remove_last_dim()-v1_pos.remove_last_dim();

        vec3 value = vec<MGLfloat,3>(
                u1[1]*u2[2] - u1[2]*u2[1],
                u1[2]*u2[0] - u1[0]*u2[2],
                u1[0]*u2[1] - u1[1]*u2[0]);
        return value.normalized();
    }

    // get Barycentric coordinates
    vec3 get_Bary(vec3 point) {
        vec4 v1_pos = this->v1.pos / this->v1.pos[3];
        vec4 v2_pos = this->v2.pos / this->v2.pos[3];
        vec4 v3_pos = this->v3.pos / this->v3.pos[3];
        
        vec3 x = v2_pos.remove_last_dim()-v1_pos.remove_last_dim();
        vec3 y = v3_pos.remove_last_dim()-v1_pos.remove_last_dim();
        vec3 z = point-v1_pos.remove_last_dim();
        MGLfloat a,b,c,area_tri; // alpha, beta, gamma and area of full triangle
        // Use of Cramer's Rule for beta and gamma since point-alpha = beta * (b-a) + gamma * (c-a)
        // From the slide 10 of "Triangles and Barycentric Coordinates"
        area_tri = x[0] * y[1] - x[1] * y[0]; // det [xy] or vectors forming the triangle
        b = z[0] * y[1] - z[1] * y[0]; // det [zy] 
        c = x[0] * z[1] - x[1] * z[0]; // det [xz]
        b /= area_tri; 
        c /= area_tri; 
        a = 1.0f - b - c;

        // First try where I got the area's as all positive
        /*area_tri = area(v1.pos.remove_last_dim(),v2.pos.remove_last_dim(),v3.pos.remove_last_dim());
        a = area(point, v2.pos.remove_last_dim(),v3.pos.remove_last_dim());
        b = area(v1.pos.remove_last_dim(), point,v3.pos.remove_last_dim());
        c = area(point, v1.pos.remove_last_dim(),v2.pos.remove_last_dim());

        a /= area_tri;
        b /= area_tri;
        c /= area_tri;*/
        //c = 1.0f - a - b;
        
        
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

// matrix structures
MGLmatrix_mode matrix_mode;
stack<mat4,vector<mat4> > projectionMatrices;
stack<mat4,vector<mat4> > modelviewMatrices;

// helper function to access current/active matrix
stack<mat4,vector<mat4> > * getCurrentStack() {
    stack<mat4,vector<mat4> > *value = NULL;
    switch(matrix_mode) {
        case MGL_MODELVIEW:
            value = &modelviewMatrices;
            break;
        case MGL_PROJECTION:
            value = &projectionMatrices;
            break;
    }
    
    // check if empty
    if (value->empty()) {
        mat4 matrix;
        matrix.make_zero();
        matrix(0,0) = 1;
        matrix(1,1) = 1;
        matrix(2,2) = 1;
        matrix(3,3) = 1;
        value->push(matrix);
    }

    return value;
}

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
MGLfloat *zbuffer = NULL;

void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    Triangle curr_triangle;
    Vertex v1, v2, v3;
    MGLfloat xmin, xmax, ymin, ymax, all;
    
    // z-buffer
    zbuffer = new MGLfloat[width*height];
    for (unsigned i = 0; i < width*height; ++i) {
        zbuffer[i] = 2.0f;
    }

    for (unsigned tri = 0; tri < tri_list.size(); tri++) {
        curr_triangle = tri_list.at(tri);
        v1 = curr_triangle.v1;
        v2 = curr_triangle.v2;
        v3 = curr_triangle.v3;

        // perspective division
        v1.pos /= v1.pos[3];
        v2.pos /= v2.pos[3];
        v3.pos /= v3.pos[3];
        
        // transform the curr_triangle's vertices to screen space in v1,v2,v3
        v1.pos[0] = (width * (v1.pos[0] + 1)) / 2;
        v1.pos[1] = (height * (v1.pos[1] + 1)) / 2;
        v2.pos[0] = (width * (v2.pos[0] + 1)) / 2;
        v2.pos[1] = (height * (v2.pos[1] + 1)) / 2;
        v3.pos[0] = (width * (v3.pos[0] + 1)) / 2;
        v3.pos[1] = (height * (v3.pos[1] + 1)) / 2;

        // bounded box of triangle
        xmin = max(0.0f,x_min(v1,v2,v3)); // shouldn't be lower than 0
        xmax = min((MGLfloat)(width-1), x_max(v1,v2,v3)); // no greater than width
        ymin = max(0.0f,y_min(v1,v2,v3)); // shouldn't be lower than 0
        ymax = min((MGLfloat)(height-1), y_max(v1,v2,v3)); // no greater than height

        // after found bounds retake the vertices from triangle
        v1 = curr_triangle.v1;
        v2 = curr_triangle.v2;
        v3 = curr_triangle.v3;

        // perspective division
        v1.pos /= v1.pos[3];
        v2.pos /= v2.pos[3];
        v3.pos /= v3.pos[3];
        
        // for each i and j in bounded box, move back to obj space and check if hits triangle
        vec3 proj_pt, bary, color; // projected point to its plane,  barycentric coords, and interpolated color
        for (MGLint i = xmin; i <= xmax; i++) {
            for (MGLint j = ymin; j <= ymax; j++) {
                // i,j --> proj_pt in obj space
                proj_pt = vec<MGLfloat,3>( (((2 *(MGLfloat)i) + 1) / width) - 1,
                        (((2 *(MGLfloat)j) + 1) / height) - 1 , 0);

                // start barycentric stuff
                bary = curr_triangle.get_Bary(proj_pt);
                
                
                proj_pt[2] = bary[0] * v1.pos[2] + bary[1] * v2.pos[2] + bary[2] * v3.pos[2];

                // perspective correct
                bary[0] /= curr_triangle.v1.pos[3];
                bary[1] /= curr_triangle.v2.pos[3];
                bary[2] /= curr_triangle.v3.pos[3];
                all = bary[0] + bary[1] + bary[2];
                bary[0] /= all;
                bary[1] /= all;
                bary[2] /= all;


                if (proj_pt[2] < -1.0f || proj_pt[2] > 1.0f || proj_pt[2] > zbuffer[width*j + i]) {
                    continue;
                }
                
                // if inside the triangle
                if (bary[0] >= 0 && bary[1] >= 0 && bary[2] >= 0) {

                    color = (bary[0] * 255.0f) * curr_triangle.v1.col + bary[1] * 255.0f * curr_triangle.v2.col + bary[2] * 255.0f * curr_triangle.v3.col;
                    data[width*j + i] = Make_Pixel(color[0], color[1], color[2]);
                    //data[width*j + i] = Make_Pixel(255,255,255);
                    zbuffer[width*j + i] = proj_pt[2];
                }
            }
        }
    }
    delete [] zbuffer;
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
    vert_list.clear();
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

    if (!modelviewMatrices.empty()){
        pos = modelviewMatrices.top() * pos;
    }
    if (!projectionMatrices.empty()) {
        pos = projectionMatrices.top() * pos;
    }
    cout << endl;

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
    
    if (!modelviewMatrices.empty()){
        pos = modelviewMatrices.top() * pos;
    }
    if (!projectionMatrices.empty()) {
        pos = projectionMatrices.top() * pos;
    }
    
    vec3 col = curr_col; 
    vert_list.push_back(Vertex(pos,col));
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    matrix_mode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
    getCurrentStack()->push(mat4(getCurrentStack()->top()));
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    getCurrentStack()->pop();
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    mat4 matrix;
    matrix.make_zero();
    for (unsigned i = 0; i < 4; i++) {
        matrix(i,i) = 1;
    }
    getCurrentStack()->top() = matrix;
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
    mat4 m1;
    for (int i = 0; i < 16; i++) {
        m1.values[i] = *(matrix+i);
    }
    getCurrentStack()->top() = m1;
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
    mat4 top = getCurrentStack()->top();
    mat4 m1;
    int m = top.rows();
    int n = top.cols();
    for (int i = 0; i < m*n; i++) {
        m1.values[i] = matrix[i];
    }
    m1 = top * m1;
    getCurrentStack()->top() = m1;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mat4 matrix;
    matrix.make_zero();
    matrix(0,0) = 1;
    matrix(1,1) = 1;
    matrix(2,2) = 1;
    matrix(3,3) = 1;
    matrix(0,3) = x;
    matrix(1,3) = y;
    matrix(2,3) = z;

    mglMultMatrix(matrix.values);
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
    // first normalze the rotation vector
    vec3 xyz = vec<MGLfloat,3>(x,y,z);
    mat4 matrix;
    MGLfloat c,s;

    xyz = xyz.normalized();
    x = xyz[0], y = xyz[1], z = xyz[2];
    angle = (M_PI * angle) / 180.0f;
    c = cos(angle); // cosine of angle
    s = sin(angle); // sine of angle
    
    // build matrix to multiply with (right to left, then down)
    matrix.make_zero();
    matrix(0,0) = x*x*(1-c) + c; 
    matrix(0,1) = x*y*(1-c) - z*s;
    matrix(0,2) = x*z*(1-c) + y*s;
    matrix(1,0) = y*x*(1-c) + z*s;
    matrix(1,1) = y*y*(1-c) + c;
    matrix(1,2) = y*z*(1-c) - x*s;
    matrix(2,0) = x*z*(1-c) - y*s;
    matrix(2,1) = y*z*(1-c) + x*s;
    matrix(2,2) = z*z*(1-c) + c;
    matrix(3,3) = 1;

    mglMultMatrix(matrix.values);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mat4 matrix;
    matrix.make_zero();
    matrix(0,0) = x;
    matrix(1,1) = y;
    matrix(2,2) = z;
    matrix(3,3) = 1.0f;
    mglMultMatrix(matrix.values);
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
    mat4 matrix;
    matrix.make_zero();
    MGLfloat rl = right - left, tb = top - bottom, fn = far - near;
    matrix(0,0) = 2.0f * near / rl;
    matrix(1,1) = 2.0f * near / tb;
    matrix(2,2) = -(far + near) / fn;
    matrix(0,2) = (right + left) / rl;
    matrix(1,2) = (top + bottom) / tb;
    matrix(3,2) = -1.0f;
    matrix(2,3) = -(2.0f * far * near) / fn;
    mglMultMatrix(matrix.values);
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
    mat4 matrix;
    matrix.make_zero();
    MGLfloat rl = right - left, tb = top - bottom, fn = far - near;
    matrix(0,0) = 2.0f / rl;
    matrix(1,1) = 2.0f / tb;
    matrix(2,2) = -2.0f / fn;
    matrix(3,3) = 1.0f;
    matrix(0,3) = -(right + left) / rl;
    matrix(1,3) = -(top + bottom) / tb;
    matrix(2,3) = -(far + near) / fn;
    
    mglMultMatrix(matrix.values);
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
