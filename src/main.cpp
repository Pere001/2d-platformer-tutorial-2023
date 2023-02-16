

#define CREATE_CONSOLE false // Set this to true to get a console you can use for debugging with DebugPrint().

#include <stdio.h>
#include <Windows.h>

#include <stdarg.h>

#include <wingdi.h>
#include <GL/gl.h>

//
// Types
//
#include <stdint.h>


//
// Define utilities
//
#define Assert(x) { if (!(x)) { *(int *)0 = 0; } }
#define AssertRange(x0, x1, x2) Assert(((x0) <= (x1)) && ((x1) <= (x2)))
#define ArrayCount(arr) (sizeof(arr) / sizeof(arr[0]))



//==============================================================================
// 
//                                  Maths
// 
//==============================================================================

#define MAX_FLOAT   340282346638528859811704183484516925440.f
#define MIN_FLOAT (-340282346638528859811704183484516925440.f)
#define PI 3.141592653589793238f
#define SQUARE(x) ((x)*(x))

#include <math.h>
inline float Cos(float x){
    float result = cos(x);
    return result;
}
inline float Sin(float x){
    float result = sin(x);
    return result;
}
inline float SquareRoot(float x){
    float result = sqrt(x);
    return result;
}

// Interpolation
inline float Lerp(float a, float b, float t){
    float result = (1.f - t)*a + b*t;
    return result;
}

// Modulo
inline float FMod(float x, float modBy){
    Assert(modBy >= 0.f);
    float result = fmod(x, modBy);
    return result;
}

// Round towards -infinity
inline float Floor(float x){
    float result = floorf(x);
    return result;
}

// Round towards +infinity
inline float Ceil(float x){
    float result = ceilf(x);
    return result;
}

inline float Round(float x){
    float result = roundf(x);
    return result;
}

// Fractional part
inline float Frac(float x){
    float result = x - (float)(int)x;
    return result;
}

// Max
inline float Max(float a, float b){
    if (a >= b)
        return a;
    return b;
}
inline int MaxInt(int a, int b){
    if (a >= b)
        return a;
    return b;
}

// Min
inline float Min(float a, float b){
    if (a <= b)
        return a;
    return b;
}
inline int MinInt(int a, int b){
    if (a <= b)
        return a;
    return b;
}

// Clamp
inline float Clamp(float value, float min, float max){
    if (value > max) return max;
    if (value < min) return min;
    return value;
}
inline float Clamp01(float value){
    return Clamp(value, 0, 1.f);
}
inline int ClampInt(int value, int min, int max){
    if (value > max) return max;
    if (value < min) return min;
    return value;
}

// Absolute value
inline float Abs(float value){
    if (value < 0)
        return -value;
    return value;
}
inline int AbsInt(int value){
    if (value < 0)
        return -value;
    return value;
}


// Return range: [-pi, pi]
inline float NormalizeAngle(float a){
    a = fmod(a, 2*PI);
    if (a < -PI)
        a += 2*PI;
    else if (a > PI)
        a -= 2*PI;
    return a;
}
// Return range. [-pi, pi]
float AngleDifference(float to, float from){
    float result = NormalizeAngle(to) - NormalizeAngle(from);
    if (result > PI)
        result -= 2*PI;
    else if (result <= -PI)
        result += 2*PI;
    return result;
}


//
// V2
//
struct v2{
    union {
        struct{
            float x;
            float y;
        };
        float asArray[2]; // TODO Get rid of this to simplify if possible.
    };
};

inline v2 V2(float x, float y){
    v2 result = {x, y};
    return result;
}
inline v2 V2(float xy){
    v2 result = {xy, xy};
    return result;
}

inline v2 operator+(v2 a, v2 b){
    v2 result = {a.x + b.x, a.y + b.y};
    return result;
}

inline v2 operator-(v2 a, v2 b){
    v2 result = {a.x - b.x, a.y - b.y};
    return result;
}

inline v2 operator-(v2 a){
    v2 result = {-a.x, -a.y};
    return result;
}

inline v2 operator/(v2 a, float scalar){
    v2 result = {a.x/scalar, a.y/scalar};
    return result;
}
inline v2 operator*(v2 a, float scalar){
    v2 result = {a.x*scalar, a.y*scalar};
    return result;
}
inline v2 operator/(float scalar, v2 a){
    v2 result = {scalar/a.x, scalar/a.y};
    return result;
}
inline v2 operator*(float scalar, v2 a){
    v2 result = {a.x*scalar, a.y*scalar};
    return result;
}

inline void operator+=(v2 &a, v2 b){
    a = a + b;
}
inline void operator*=(v2 &a, float scalar){
    a = a * scalar;
}
inline void operator-=(v2 &a, v2 b){
    a = a - b;
}
inline void operator/=(v2 &a, float scalar){
    a = a / scalar;
}
inline bool operator==(v2 a, v2 b){
    bool result = (a.x == b.x) && (a.y == b.y);
    return result;
}
inline bool operator!=(v2 a, v2 b){
    bool result = (a.x != b.x) || (a.y != b.y);
    return result;
}

inline float Dot(v2 a, v2 b){
    float result = a.x*b.x + a.y*b.y;
    return result;
}

inline float Cross(v2 a, v2 b){
    float result = a.x*b.y - a.y*b.x;
    return result;
}

inline v2 Hadamard(v2 a, v2 b){
    v2 result = {a.x*b.x, a.y*b.y};
    return result;
}

inline float Length(v2 a){
    float result = SquareRoot(a.x*a.x + a.y*a.y);
    return result;
}
inline float LengthSqr(v2 a){
    float result = a.x*a.x + a.y*a.y;
    return result;
}

inline float AngleOf(v2 a){
    float result = 0;
    if (a.x || a.y)
        result = atan2(a.y, a.x);
    return result;
}

inline v2 V2FromLengthDir(float length, float direction){
    v2 result = {Cos(direction)*length, Sin(direction)*length};
    return result;
}
// Counterclockwise
inline v2 Rotate90Degrees(v2 p){
    v2 result = {-p.y, p.x};
    return result;
}
// Counterclockwise
inline v2 RotateMinus90Degrees(v2 p){
    v2 result = {p.y, -p.x};
    return result;
}

inline v2 Normalize(v2 a){
    v2 result = V2FromLengthDir(1.f, AngleOf(a));
    return result;
}

inline bool PointInCircle(v2 p, v2 c, float r){
    float distanceSqr = LengthSqr(p - c);
    if (distanceSqr <= r*r)
        return true;
    return false;
}
// - Tangent doesn't count as inside.
inline bool PointInRectangle(v2 point, v2 rectMin, v2 rectMax){
    bool result = (point.x > rectMin.x & point.x < rectMax.x & point.y > rectMin.y & point.y < rectMax.y);
    return result;
}

inline v2 FloorV2(v2 a){
    v2 result = {Floor(a.x), Floor(a.y)};
    return result;
}
inline v2 MinV2(v2 a, v2 b){
    v2 result = {Min(a.x, b.x), Min(a.y, b.y)};
    return result;
}
inline v2 MaxV2(v2 a, v2 b){
    v2 result = {Max(a.x, b.x), Max(a.y, b.y)};
    return result;
}

bool PointsAreCW(v2 *points, int num){
    float angleSum = 0.f;
    v2 prevPoint = points[num - 1];
    float prevAngle = AngleOf(prevPoint - points[num - 2]);
    for(int i = 0; i < num; i++){
        float newAngle = AngleOf(points[i] - prevPoint);
        angleSum += AngleDifference(newAngle, prevAngle);
        prevAngle = newAngle;
        prevPoint = points[i];
    }
    return (angleSum < 0);
}



//==============================================================================
// 
//                            Some Platform Stuff
// 
//==============================================================================

struct button_state {
    bool isDown;
    int transitionCount;
};

bool ButtonWentDown(button_state *b){
    if (b->isDown && (b->transitionCount % 2)){
        return true;
    }
    return false;
}
bool ButtonWentUp(button_state *b){
    if (!b->isDown && (b->transitionCount % 2)){
        return true;
    }
    return false;
}

struct keyboard_input{
    union{
        button_state asArray[48];
        struct{
            button_state letters[26];
            button_state numbers[10];
            button_state escape;
            button_state enter;
            button_state space;
            button_state shift;
            button_state control;
            button_state backspace;
            button_state alt;
            button_state tab;
            button_state arrowLeft;
            button_state arrowRight;
            button_state arrowUp;
            button_state arrowDown;
        };
    };
};
struct input_state {
    keyboard_input keyboard;
    button_state mouseButtons[5];
    v2 mousePos;
    v2 windowDim;
};

void UpdateButtonState(button_state *b, bool newIsDown){
    if (!b->isDown != !newIsDown){
        b->isDown = newIsDown;
        b->transitionCount++;
    }
}

static HANDLE globalStdHandle = {};
static LARGE_INTEGER globalPerformanceFrequency;
static bool globalRunning;
static input_state globalInput = {};
static HCURSOR globalCursor;
static char globalExePath[1024] = {};

void DebugPrint(char *str){
    //OutputDebugStringA(str);
    WriteFile(globalStdHandle, str, (DWORD)strlen(str), 0, 0);
}
void DebugPrintf(char *format, ...){
    va_list argptr;
    va_start(argptr, format);

    char localStr[1024];
    vsprintf_s(localStr, sizeof(localStr), format, argptr);

    va_end(argptr);

    localStr[1023] = 0; // null-terminate
    DebugPrint(localStr);
}



//==============================================================================
//
//                                   Game State
//
//==============================================================================

struct wall{
    v2 p[3];
    v2 normals[3];
};
// Create wall from 3 points.
wall Wall(v2 p0, v2 p1, v2 p2){
    wall w;
    w.p[0] = p0;
    w.p[1] = p1;
    w.p[2] = p2;

    if (Cross(w.p[1] - w.p[0], w.p[2] - w.p[1]) < 0){ // Points are CW
        // Make them CCW
        v2 temp = w.p[1];
        w.p[1] = w.p[2];
        w.p[2] = temp;
    }

    // Compute normals
    v2 edgeDir01 = Normalize(w.p[1] - w.p[0]);
    v2 edgeDir12 = Normalize(w.p[2] - w.p[1]);
    v2 edgeDir20 = Normalize(w.p[0] - w.p[2]);
    w.normals[0] = V2(edgeDir01.y, -edgeDir01.x);
    w.normals[1] = V2(edgeDir12.y, -edgeDir12.x);
    w.normals[2] = V2(edgeDir20.y, -edgeDir20.x);

    return w;
}
struct game_state{
    v2 cameraPos;
    float cameraScale;

    v2 rectPos;
    v2 rectDim;
    v2 rectSpeed;

    v2 circlePos;
    float circleRadius;
    v2 circleSpeed;

    int numWalls;
    wall walls[100];
};
static game_state gameState = {};



//==============================================================================
// 
//                           Colllision Detection
// 
//==============================================================================

bool PointWallCollision(v2 point, wall *w){
    for(int i = 0; i < 3; i++){
        float proj = Dot(w->normals[i], point - w->p[i]); // Center projected onto normal
        if (proj > 0) // Too far outwards from the edge
            return false;
    }
    return true;
}
bool PointWallCollision(v2 point){
    for(int i = 0; i < gameState.numWalls; i++){
        if (PointWallCollision(point, &gameState.walls[i]))
            return true;
    }
    return false;
}

bool RectangleWallCollision(v2 rMin, v2 rMax, wall *w){
    // Check if rectangle and triangle are separated by an axis-aligned line.
    // (The bitwise '&' operator here does the same as the boolean '&&', but it's faster because it avoids some branches)
    if (   (w->p[0].x >= rMax.x & w->p[1].x >= rMax.x & w->p[2].x >= rMax.x)
        || (w->p[0].y >= rMax.y & w->p[1].y >= rMax.y & w->p[2].y >= rMax.y)
        || (w->p[0].x <= rMin.x & w->p[1].x <= rMin.x & w->p[2].x <= rMin.x)
        || (w->p[0].y <= rMin.y & w->p[1].y <= rMin.y & w->p[2].y <= rMin.y))
        return false;

    // Check if rectangle and triangle are separated by one of triangle's edges.
    for(int i = 0; i < 3; i++){
        v2 p = w->p[i];
        v2 n = w->normals[i];

        float rectProjected[4];
        rectProjected[0] = Dot(rMin               - p, n);
        rectProjected[1] = Dot(V2(rMin.x, rMax.y) - p, n);
        rectProjected[2] = Dot(rMax               - p, n);
        rectProjected[3] = Dot(V2(rMax.x, rMin.y) - p, n);

        float rectProjectedMin = Min(Min(Min(rectProjected[0], rectProjected[1]), rectProjected[2]), rectProjected[3]);

        if (rectProjectedMin >= 0){
            return false;
        }
    }
    return true; // No separating axis found
}
bool RectangleWallCollision(v2 rMin, v2 rMax){
    for(int i = 0; i < gameState.numWalls; i++){
        if (RectangleWallCollision(rMin, rMax, &gameState.walls[i]))
            return true;
    }
    return false;
}

bool CircleWallCollision(v2 c, float r, wall *w){
    // This loop checks the collision of the center with the triangle and the edge rectangles.
    bool centerOnTriangle = true;
    for(int i = 0; i < 3; i++){
        float proj = Dot(w->normals[i], c - w->p[i]); // Center projected onto normal
        
        if (proj > r) // (Optional early-out: circle is too far outwards to collide)
            return false;

        if (proj < 0) // The center is inward from this edge: it might be inside the triangle or not,
            continue; // and inside other edge rects or not, but this will be found by the other iterations.

        // Check if center is inside the edge rects.
        v2 edgeDir = Rotate90Degrees(w->normals[i]);
        float edgeProj0 = Dot(edgeDir, c - w->p[i]); // Center projected onto edge direction, relative to current vertex.
        float edgeProj1 = Dot(edgeDir, c - w->p[(i + 1) % 3]); // Center projected onto edge direction, relative to next vertex.
        if (edgeProj0 > 0 && edgeProj1 < 0){
            // The projected center falls between the two vertices of the edge: the center is inside the edge rect.
            return true;
        }
        
        centerOnTriangle = false;
    }

    if (!centerOnTriangle){
        // Center wasn't inside the triangle or the edge rects
        // Check if a vertex is inside the circle.
        for(int i = 0; i < 3; i++){
            float distanceSqr = LengthSqr(w->p[i] - c);
            if (distanceSqr < r*r)
                return true;
        }
    }
    return centerOnTriangle;
}
bool CircleWallCollision(v2 c, float r){
    for(int i = 0; i < gameState.numWalls; i++){
        if (CircleWallCollision(c, r, &gameState.walls[i]))
            return true;
    }
    return false;
}



//==============================================================================
// 
//                                   Movement
// 
//==============================================================================

v2 MoveCircle(v2 pos, float r, v2 speed, bool *outCollided, v2 *outCollisionNormal){
    if (speed == V2(0)){
        *outCollided = false;
        return pos;
    }
    
    wall *walls = gameState.walls;
    int numWalls = gameState.numWalls;

    //
    // Move
    //
    bool collided = false;

    v2 newPos = pos;
    float numerator = 1.f;
    float denominator = 1.f;

    v2 lastCollisionPos = {};
    int firstCollidingWallIndex = numWalls; // (This default value will cause the collision normal code to be skipped.)
    
    bool alreadyShifted = false; // We can only shift once
    v2 permanentShift = {0, 0}; // After having shifted, all subsequent steps will apply the same shift, stored here.

    float speedLength = Length(speed);

    float stepSize = 1.f;
    float speedFactor = 1.f;

    for(int steps = 0; steps < 15; steps++){
        v2 p = pos + speedFactor*speed + permanentShift;

        bool collision = false;
        for(int i = 0; i < numWalls; i++){ // Check collision against all walls
            wall *wll = &walls[i];
            if (CircleWallCollision(p, r, wll)){
                if (alreadyShifted){ // Accept the collision
                    collision = true;
                    lastCollisionPos = p;
                    firstCollidingWallIndex = i;
                }else{ // Try to shift
                    // "Shifting" means we'll check again for collision after applying a tiny offset in both directions
                    // perpendicular to the speed, in order to facilitate smooth movement tangent to a wall.

                    float shiftDistance = .1f;
                    v2 shiftVector = shiftDistance*Normalize(Rotate90Degrees(speed));
                    v2 shifts[2] = {shiftVector, -shiftVector};
                    bool shiftFailed[2] = {false, false};

                    // We store p + shift here to avoid it having different results when calculated in different places
                    // because of optimizer or something, which could potentially put us inside a wall.
                    v2 shiftedP[2] = {p + shifts[0], p + shifts[1]};
                    
                    for(int j = 0; j < numWalls; j++){
                        wall *w = &walls[j];
                        for(int k = 0; k < 2; k++){
                            if (!shiftFailed[k] && CircleWallCollision(shiftedP[k], r, w)){
                                shiftFailed[k] = true;
                            }
                        }    
                        if (shiftFailed[0] & shiftFailed[1])
                            break;
                    }
                    
                    if (shiftFailed[0] & shiftFailed[1]){ // Both shifts were unsuccessful
                        collision = true;
                        lastCollisionPos = p;
                        firstCollidingWallIndex = i;
                    }else{
                        // One of the shifts is a free position.
                        for(int j = 0; j < 2; j++){
                            if (!shiftFailed[j]){
                                p = shiftedP[j]; // p will be assigned to newPos.
                                permanentShift = shifts[j];
                                alreadyShifted = true;
                                break;
                            }
                        }
                    }
                    break; // Stop iterating walls for this step (because we either found a collision or already iterated all walls and found a free shifted position).
                }
            }
        }

        stepSize *= .5f;
        if (collision){
            speedFactor -= stepSize;
            collided = true;
        }else{
            speedFactor += stepSize;
            newPos = p;
            if (steps == 0)
                break; // No collisions found on the first step.
        }

        if (stepSize*speedLength <= .05f) // Limit step length. 
            break;
    }
    

    //
    // Calculate normal
    //
    v2 bestNormal = -Normalize(speed); // Default value in case we don't find any.
    float bestTime = MAX_FLOAT; // Lowest time until collision, to find the edge that would collide first.
    for(int wallIndex = firstCollidingWallIndex; wallIndex < numWalls; wallIndex++){
        wall *w = &walls[wallIndex];
        if (!CircleWallCollision(lastCollisionPos, r, w))
            continue;

        for(int i = 0; i < 3; i++){
            if (Dot(w->normals[i], speed) >= 0)
                continue; // Ignore edges facing away

            v2 t0 = w->p[i];
            v2 t1 = w->p[(i + 1) % 3];
        
            // Let's find the closest point on the edge.
            v2 edgeDir = Rotate90Degrees(w->normals[i]);
            float cProjected = Dot(edgeDir, newPos - t0); // Center projected onto the edge.
            float t1Projected = Dot(edgeDir, t1 - t0); // Next vertex projected onto the edge.
            v2 collisionPoint = t0 + edgeDir*Clamp(cProjected, 0, t1Projected); // Closest point on the edge

            v2 n = Normalize(newPos - collisionPoint); // Collision normal of this edge.
            float distance = Dot(newPos - collisionPoint, n) - r; // (We use Dot() to avoid the squareroot of Length()).
            float nDotSpeed = Dot(n, speed);
            if (!nDotSpeed)
                continue; // Avoid division by 0.
            float time = distance/-nDotSpeed; // Time to hit wall
            if (time < bestTime && Dot(n, speed) <= 0){ // Ignore normals that face away from speed (rare but technically possible because of shifts)
                bestTime = time;
                bestNormal = n;
            }
        }
    }

    *outCollided = collided;
    *outCollisionNormal = bestNormal;
    return newPos;
}


// Collision box is centered at pos.
v2 MoveRectangle(v2 pos, v2 halfDim, v2 speed, bool *outCollided, v2 *outCollisionNormal){
    if (speed == V2(0)){
        *outCollided = false;
        return pos;
    }

    wall *walls = gameState.walls;
    int numWalls = gameState.numWalls;

    //
    // Move
    //
    bool collided = false;

    v2 newPos = pos;

    v2 lastCollisionPos = {};
    int firstCollidingWallIndex = numWalls; // (This default value will cause the collision normal code to be skipped.)
    
    bool alreadyShifted = false; // We can only shift once
    v2 permanentShift = {0, 0}; // After having shifted, all subsequent steps will apply the same shift, stored here.

    float speedLength = Length(speed);

    float stepSize = 1.f;
    float speedFactor = 1.f;

    for(int steps = 0; steps < 15; steps++){
        v2 p = pos + speedFactor*speed + permanentShift;

        bool collision = false;
        for(int i = 0; i < numWalls; i++){ // Check collision against all walls
            wall *wll = &walls[i];
            if (RectangleWallCollision(p - halfDim, p + halfDim, wll)){
                if (alreadyShifted){ // Accept the collision
                    collision = true;
                    firstCollidingWallIndex = i;
                }else{ // Try to shift
                    // "Shifting" means we'll check again for collision after applying a tiny offset in both directions
                    // perpendicular to the speed, in order to facilitate smooth movement tangent to a wall.
                    
                    float shiftDistance = .1f;
                    v2 shiftVector = shiftDistance*Normalize(Rotate90Degrees(speed));
                    v2 shifts[2] = {shiftVector, -shiftVector};
                    bool shiftFailed[2] = {false, false};
                    
                    // We store p + shift here to avoid it having different results when calculated in different places
                    // because of optimizer or something, which could potentially put us inside a wall.
                    v2 shiftedP[2] = {p + shifts[0], p + shifts[1]};
                    
                    for(int j = 0; j < numWalls; j++){
                        wall *w = &walls[j];
                        for(int k = 0; k < 2; k++){
                            if (!shiftFailed[k] && RectangleWallCollision(shiftedP[k] - halfDim, shiftedP[k] + halfDim, w)){
                                shiftFailed[k] = true;
                            }
                        }    
                        if (shiftFailed[0] & shiftFailed[1])
                            break;
                    }
                    
                    if (shiftFailed[0] & shiftFailed[1]){ // Both shifts were unsuccessful
                        collision = true;
                        firstCollidingWallIndex = i;
                    }else{
                        // One of the shifts is a free position.
                        for(int j = 0; j < 2; j++){
                            if (!shiftFailed[j]){
                                p = shiftedP[j]; // p will be assigned to newPos.
                                permanentShift = shifts[j];
                                alreadyShifted = true;
                                break;
                            }
                        }
                    }
                    break; // Stop iterating walls for this step (because we either found a collision or already iterated all walls and found a free shifted position).
                }
            }
        }

        stepSize *= .5f;
        if (collision){
            speedFactor -= stepSize;
            collided = true;
            lastCollisionPos = p;
        }else{
            speedFactor += stepSize;
            newPos = p;
            if (steps == 0)
                break; // No collisions found on the first step.
        }

        if (stepSize*speedLength <= .05f) // Limit step length. 
            break;
    }
    
    //
    // Calculate normal
    //
    
    v2 newR0 = newPos - halfDim;
    v2 newR1 = newPos + halfDim;
    v2 newRectPoints[4] = {newR0, V2(newR0.x, newR1.y), newR1, V2(newR1.x, newR0.y)};

    // We find the two edges (one horizontal and one vertical) that face the direction of the speed,
    // and these arrays hold their normal and a point of the edge.
    v2 rectEdgePoints[2] = { {(speed.x < 0 ? newR0.x : newR1.x), newPos.y},
                             {newPos.x, (speed.y < 0 ? newR0.y : newR1.y)} };
    v2 rectEdgeNormals[2] = { {(speed.x < 0 ? -1.f : 1.f), 0},
                              {0, (speed.y < 0 ? -1.f : 1.f)}};

    v2 bestNormal = -Normalize(speed); // Final collision normal (Default value in case we don't find any)
    float bestTime = MAX_FLOAT; // Lowest time until collision, to find the edge that would collide first.

    for(int wallIndex = 0; wallIndex < numWalls; wallIndex++){
        wall *w = &walls[wallIndex];
        if (!RectangleWallCollision(lastCollisionPos - halfDim, lastCollisionPos + halfDim, w))
            continue;
    
        v2 localNormal = bestNormal; // Best normal of the current wall. (default value "just in case")
        float localBestTime = MIN_FLOAT; // Highest time to hit an edge

        // Find the edge with the highest time until impact
    
        // Project rectangle onto triangle normals
        for(int i = 0; i < 3; i++){
            v2 n = w->normals[i];
            v2 p = w->p[i];
            float nDotSpeed = Dot(n, speed);
            if (nDotSpeed < 0){ // Ignore edges facing away.
                float d = Min(Min(Min(Dot(newRectPoints[0] - p, n), // Projected distance to rectangle
                                      Dot(newRectPoints[1] - p, n)),
                                      Dot(newRectPoints[2] - p, n)),
                                      Dot(newRectPoints[3] - p, n));
                float t = d/-nDotSpeed; // Time to hit edge
                if (t > localBestTime){
                    localBestTime = t;
                    localNormal = n;
                }
            }
        }
        
        // Project triangle onto rectangle normals (we only need to check the two edges facing the direction of the speed)
        for(int i = 0; i < 2; i++){
            v2 n = rectEdgeNormals[i];
            float nDotSpeed = Dot(n, speed);
            if (nDotSpeed){ // Gota do this to avoid potential division by 0.
                float d = Min(Min(Dot(w->p[0] - rectEdgePoints[i], n), // Projected distance to triangle
                                  Dot(w->p[1] - rectEdgePoints[i], n)),
                                  Dot(w->p[2] - rectEdgePoints[i], n));
                float t = d/nDotSpeed; // Time to hit edge
                if (t > localBestTime){
                    localBestTime = t;
                    localNormal = -n;
                }
            }
        }

        // Compare the best normal in this wall with the best of previous walls.
        float epsilon = 0.000001f;
        if (localBestTime < bestTime ||
            (localBestTime - bestTime < epsilon && Dot(localNormal, speed) < Dot(bestNormal, speed))) // If there's a tie we'll take the normal most aligned with the speed.
        {
            bestNormal = localNormal;
            bestTime = localBestTime;
        }
    }


    *outCollided = collided;
    *outCollisionNormal = bestNormal;
    return newPos;
}



//==============================================================================
// 
//                                 Platform Stuff
// 
//==============================================================================

LARGE_INTEGER GetCurrentTimeCounter(){
    LARGE_INTEGER result = {};
    QueryPerformanceCounter(&result);
    return result;
}

float GetSecondsElapsed(LARGE_INTEGER t0, LARGE_INTEGER t1){
    float result = (float)(t1.QuadPart - t0.QuadPart) / (float)globalPerformanceFrequency.QuadPart;
    return result;
}

v2 GetWindowDimension(HWND window){
    RECT rect = {};
    GetClientRect(window, &rect);

    v2 result = {(float)(rect.right - rect.left), (float)(-rect.top + rect.bottom)};
    return result;
}

LRESULT CALLBACK WindowProc(HWND window, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_DESTROY:
    {
        PostQuitMessage(0);
        globalRunning = false;
        return 0;
    }

    }

    return DefWindowProc(window, msg, wParam, lParam);
}



//==============================================================================
// 
//                                    main()
// 
//==============================================================================
int WINAPI WinMain(HINSTANCE instance, HINSTANCE prevInstance, PSTR cmdLine, int cmdShow)
{
    QueryPerformanceFrequency(&globalPerformanceFrequency);

    globalCursor = LoadCursor(0, IDC_ARROW);

    // Load exe path
    GetModuleFileNameA(0, globalExePath, ArrayCount(globalExePath));
    {
        int len = strlen(globalExePath);
        for(int i = len-1; i >= 0; i--){
            if (globalExePath[i] == '\\')
                break;
            globalExePath[i] = '\0';
        }
    }
    //
    // Create window
    //
    WNDCLASS windowClass = {};
    windowClass.lpfnWndProc   = WindowProc;
    windowClass.hInstance     = instance;
    windowClass.lpszClassName = "main window class";
    windowClass.style = CS_HREDRAW | CS_VREDRAW;
    RegisterClass(&windowClass);

    v2 initialWindowSize = {1220.f, 840.f};
    RECT windowFrameSize = {0, 0, (int)initialWindowSize.x, (int)initialWindowSize.y};
    AdjustWindowRect(&windowFrameSize, WS_OVERLAPPEDWINDOW|WS_VISIBLE, false); 
    HWND window = CreateWindowEx(0, windowClass.lpszClassName, "Physics Article Demo",
                                 WS_OVERLAPPEDWINDOW|WS_VISIBLE,
                                 CW_USEDEFAULT, CW_USEDEFAULT,
                                 windowFrameSize.right, windowFrameSize.bottom,
                                 //(int)initialWindowSize.x, (int)initialWindowSize.y,
                                 0, 0, instance, 0);
    if (!window)
        return 0;
    //ShowWindow(window, cmdShow);
    HCURSOR cursor = LoadCursor(NULL, IDC_ARROW);
    SetCursor(cursor);

    

    //
    // Create console
    //
    if (CREATE_CONSOLE){
        if(AttachConsole((DWORD)-1) == 0){ // wasn't launched from console
            AllocConsole(); // alloc your own instead
        }
    }
    globalStdHandle = GetStdHandle(STD_OUTPUT_HANDLE);


    globalInput.windowDim = GetWindowDimension(window);

    HDC dc = GetDC(window);

    //
    // OpenGl
    //
    PIXELFORMATDESCRIPTOR pfd = {};
    pfd.nSize        = sizeof(pfd);
    pfd.nVersion     = 1;
    pfd.dwFlags      = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL;
    pfd.iPixelType   = PFD_TYPE_RGBA;
    pfd.cColorBits   = 32;
    int pf = ChoosePixelFormat(dc, &pfd);

    if (!pf)
        return 0;
    
    if (!SetPixelFormat(dc, pf, &pfd))
        return 0;

    DescribePixelFormat(dc, pf, sizeof(pfd), &pfd);

    // Device context
    HGLRC rc = wglCreateContext(dc);
    wglMakeCurrent(dc, rc);

    int timeInFrames = 0;
    
    //
    // Init Game State
    //
    gameState.cameraPos = V2(0, 0);
    gameState.cameraScale = 1.f;

    gameState.rectPos = V2(0, 0);
    gameState.rectDim = V2(40, 60);
    gameState.circlePos = V2(98, 0);
    gameState.circleRadius = 50.f;

    gameState.numWalls = 20;
    gameState.walls[0] = Wall(V2(-250.f, -50.f), V2(-250.f, -150.f), V2(150.f, -50.f)); // Center Floor
    gameState.walls[1] = Wall(V2(150.f, -50.f), V2(-250.f, -150.f), V2(150.f, -150.f));

    gameState.walls[2] = Wall(V2(-250.f, -150.f), V2(-250.f, 50.f), V2(-350.f, 50.f)); // Left Floor
    gameState.walls[3] = Wall(V2(-250.f, -150.f), V2(-350.f, 50.f), V2(-350.f, -150.f));

    gameState.walls[4] = Wall(V2(150, -150), V2(150, -50), V2(550, -50)); // Right Floor
    gameState.walls[5] = Wall(V2(150, -150), V2(550, -50), V2(550, -150));
    
    gameState.walls[6] = Wall(V2(-350.f, 50.f), V2(-550.f, 150.f), V2(-550.f, 50.f)); // Left Slope
    gameState.walls[7] = Wall(V2(-350.f, 50.f), V2(-550.f, 50.f), V2(-350.f, -150.f));
    gameState.walls[8] = Wall(V2(-115, -43), V2(-248, 33), V2(-250, -50));
    
    gameState.walls[9] = Wall(V2(150, -50), V2(250, 50), V2(250, -50)); // Right Slopes
    gameState.walls[10] = Wall(V2(250, 50), V2(350, -50), V2(250, -50));
    gameState.walls[11] = Wall(V2(350, -50), V2(550, 150), V2(550, -50));
    gameState.walls[12] = Wall(V2(321, -60), V2(386, -53), V2(351, -20)); // Intersecting wall

    gameState.walls[13] = Wall(V2(407, 117), V2(429, 270), V2(284, 150)); // Random flying walls
    gameState.walls[14] = Wall(V2(-120, 220), V2(53, 63), V2(133, 265));
    gameState.walls[15] = Wall(V2(-384, 193), V2(-350, 250), V2(-400, 250));

    gameState.walls[16] = Wall(V2(550, 150), V2(600, 150), V2(550, 500)); // Right edge barrier
    gameState.walls[17] = Wall(V2(600, 150), V2(600, 500), V2(550, 500));
    gameState.walls[18] = Wall(V2(-600, 150), V2(-550, 150), V2(-600, 500)); // Left edge barrier
    gameState.walls[19] = Wall(V2(-550, 150), V2(-550, 500), V2(-600, 500));



    LARGE_INTEGER lastFrameTime = GetCurrentTimeCounter();
    globalRunning = true;
    while(globalRunning){

        //
        // Message Loop
        //
        MSG msg = { };
        while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE) > 0) {
            switch(msg.message){
            case WM_KEYUP:
            case WM_KEYDOWN:
            {
                //
                // Keyboard Input
                //
                bool wentDown = (msg.message == WM_KEYDOWN);
                auto k = &globalInput.keyboard;
                if (msg.wParam == VK_ESCAPE){
                    UpdateButtonState(&k->escape, wentDown);
                }else if (msg.wParam == VK_RETURN){
                    UpdateButtonState(&k->enter, wentDown);
                }else if (msg.wParam == VK_SPACE){
                    UpdateButtonState(&k->space, wentDown);
                }else if (msg.wParam == VK_SHIFT){
                    UpdateButtonState(&k->shift, wentDown);
                }else if (msg.wParam == VK_CONTROL){
                    UpdateButtonState(&k->control, wentDown);
                }else if (msg.wParam == VK_BACK){
                    UpdateButtonState(&k->backspace, wentDown);
                }else if (msg.wParam == VK_MENU){
                    UpdateButtonState(&k->alt, wentDown);
                }else if (msg.wParam == VK_TAB){
                    UpdateButtonState(&k->tab, wentDown);
                }else if (msg.wParam == VK_LEFT){
                    UpdateButtonState(&k->arrowLeft, wentDown);
                }else if (msg.wParam == VK_RIGHT){
                    UpdateButtonState(&k->arrowRight, wentDown);
                }else if (msg.wParam == VK_UP){
                    UpdateButtonState(&k->arrowUp, wentDown);
                }else if (msg.wParam == VK_DOWN){
                    UpdateButtonState(&k->arrowDown, wentDown);
                }else if (msg.wParam >= 'A' && msg.wParam <= 'Z'){
                    UpdateButtonState(&k->letters[msg.wParam - 'A'], wentDown);
                }else if (msg.wParam >= '0' && msg.wParam <= '9'){
                    UpdateButtonState(&k->numbers[msg.wParam - '0'], wentDown);
                }

            } break;
            
            default:
            {
                TranslateMessage(&msg);
                DispatchMessage(&msg);
            }
            }
        }


        //
        // Mouse Input
        //
        globalInput.windowDim = GetWindowDimension(window);

        POINT mousePoint;
        GetCursorPos(&mousePoint);
        ScreenToClient(window, &mousePoint);
        globalInput.mousePos.x = (float)(int)mousePoint.x;
        globalInput.mousePos.y = (float)(int)(globalInput.windowDim.y - mousePoint.y);

        bool hasFocus = (GetFocus() == window);
        if (PointInRectangle(globalInput.mousePos, V2(-.1f), globalInput.windowDim) && hasFocus){
            UpdateButtonState(&globalInput.mouseButtons[0], GetKeyState(VK_LBUTTON)  & (1 << 15));
            UpdateButtonState(&globalInput.mouseButtons[1], GetKeyState(VK_MBUTTON)  & (1 << 15));
            UpdateButtonState(&globalInput.mouseButtons[2], GetKeyState(VK_RBUTTON)  & (1 << 15));
            UpdateButtonState(&globalInput.mouseButtons[3], GetKeyState(VK_XBUTTON1) & (1 << 15));
            UpdateButtonState(&globalInput.mouseButtons[4], GetKeyState(VK_XBUTTON2) & (1 << 15));
        }

        //
        // Game Code
        //
        v2 windowDim = globalInput.windowDim;
        v2 viewPos = gameState.cameraPos - (windowDim/gameState.cameraScale)/2 + V2(0, 140);
        v2 mouseWorldPos = viewPos + globalInput.mousePos/gameState.cameraScale;


        // Left click: Place walls
        static int placingWallPointIndex = 0;
        if (!globalInput.keyboard.shift.isDown){
            if (gameState.numWalls < ArrayCount(gameState.walls)){
                wall *w = &gameState.walls[gameState.numWalls];
                v2 pointPos = mouseWorldPos;
                if (globalInput.keyboard.control.isDown){ // Snap to grid
                    float gridSize = 50;
                    pointPos = FloorV2((mouseWorldPos + V2(gridSize/2))/gridSize)*gridSize;
                }
                w->p[placingWallPointIndex] = pointPos;
                if (ButtonWentDown(&globalInput.mouseButtons[0])){
                    placingWallPointIndex = (placingWallPointIndex + 1) % 3;
                    if (placingWallPointIndex == 0){ // Made a new wall
                        // Set normals
                        gameState.walls[gameState.numWalls] = Wall(gameState.walls[gameState.numWalls].p[0], gameState.walls[gameState.numWalls].p[1], gameState.walls[gameState.numWalls].p[2]);
                        DebugPrintf("Placed a wall {(%.0f, %.0f), (%.0f, %.0f), (%.0f, %.0f)}\n",
                                    gameState.walls[gameState.numWalls].p[0].x, gameState.walls[gameState.numWalls].p[0].y,
                                    gameState.walls[gameState.numWalls].p[1].x, gameState.walls[gameState.numWalls].p[1].y,
                                    gameState.walls[gameState.numWalls].p[2].x, gameState.walls[gameState.numWalls].p[2].y);
                        gameState.numWalls++;
                    }
                    if (gameState.numWalls < ArrayCount(gameState.walls))
                        gameState.walls[gameState.numWalls].p[placingWallPointIndex] = pointPos;
                }
            }
        }
        
        // Right click: Remove walls.
        if (ButtonWentDown(&globalInput.mouseButtons[2])){
            if (placingWallPointIndex){
                placingWallPointIndex = 0; // Cancel current half-baked wall.
            }else{ // Remove walls.
                for(int i = 0; i < gameState.numWalls; i++){
                    if (PointWallCollision(mouseWorldPos, &gameState.walls[i])){
                        gameState.numWalls--;
                        if (i < gameState.numWalls){ // Fill hole in array
                            memmove(&gameState.walls[i], &gameState.walls[i + 1], sizeof(wall)*(gameState.numWalls - i));
                        }
                    }
                }
            }
        }

        // 
        // Player Movement: Circle
        //
        {
            // Read input and set speed.
            bool grounded;
            v2 groundNormal;
            // We just call this to get the grounded state
            MoveCircle(gameState.circlePos, gameState.circleRadius, V2(0, -1.f), &grounded, &groundNormal);
            if (!grounded)
                groundNormal = V2(0, 1.f);
            
            // Read input and set speed.
            // Move in direction of ground to avoid "stepping" down slopes.
            float ax = 0; // "Horizontal" acceleration
            v2 groundDir = RotateMinus90Degrees(groundNormal);
            if (globalInput.keyboard.letters['D' - 'A'].isDown){
                ax = .7f; // Accelerate right
            }else if (globalInput.keyboard.letters['A' - 'A'].isDown){
                ax = -.7f; // Accelerate left
            }else{
                ax = Dot(gameState.circleSpeed, groundDir)*-.1f; // Decelerate.
            }
            gameState.circleSpeed += ax*groundDir;
            // Limit speed
            float maxSpeed = 22.f;
            if (Length(gameState.circleSpeed) > maxSpeed)
                gameState.circleSpeed *= maxSpeed/Length(gameState.circleSpeed);


            float minSpeedY = -30.f;
            if (grounded){
                if (globalInput.keyboard.letters['W' - 'A'].isDown){
                    gameState.circleSpeed.y = 10.f; // Jump
                }else{
                    float diff = Abs(AngleDifference(PI/2, AngleOf(groundNormal)));
                    if (diff > PI/4 && gameState.circleSpeed.y > minSpeedY)
                        gameState.circleSpeed.y -= .3f;
                }
            }else if (gameState.circleSpeed.y > minSpeedY){
                gameState.circleSpeed.y -= .4f; // Gravity
            }

            float speedLength = Length(gameState.circleSpeed);
            float toMove = speedLength;
            for(int i = 0; i < 6; i++){ // Artificial iteration limit
                bool collided;
                v2 collisionNormal;
                v2 oldPos = gameState.circlePos;
                gameState.circlePos = MoveCircle(gameState.circlePos, gameState.circleRadius, V2FromLengthDir(toMove, AngleOf(gameState.circleSpeed)), &collided, &collisionNormal);
                toMove -= Length(gameState.circlePos - oldPos);
                if (!collided)
                    break;
                v2 effectiveSpeed = V2FromLengthDir(speedLength, AngleOf(gameState.circleSpeed));
                if (Dot(effectiveSpeed, collisionNormal) < -9.f) {
                    // Hit a wall hard: bounce!
                    gameState.circleSpeed = (gameState.circleSpeed - 2*Dot(gameState.circleSpeed, collisionNormal)*collisionNormal)*.4f;
                }else if (Abs(AngleDifference(AngleOf(collisionNormal), PI/2)) > PI*.2f){
                    // There's considerable slope: slide!
                    gameState.circleSpeed = gameState.circleSpeed - Dot(gameState.circleSpeed, collisionNormal)*collisionNormal;
                }else{
                    gameState.circleSpeed = gameState.circleSpeed - Dot(gameState.circleSpeed, collisionNormal)*collisionNormal;
                    break;
                }
                if (toMove < .1f)
                    break;
            }
        }


        // 
        // Player Movement: Rectangle
        //
        {
            v2 oldPos = gameState.rectPos;

            bool grounded;
            v2 groundNormal;
            // We just call this to get the grounded state
            MoveRectangle(gameState.rectPos, gameState.rectDim/2, V2(0, -1.f), &grounded, &groundNormal);
            if (!grounded){
                groundNormal = V2(0, 1.f);
            }

            // Read input and set speed.
            // Move in direction of ground to avoid "stepping" down slopes.
            float ax = 0; // "Horizontal" acceleration
            v2 groundDir = RotateMinus90Degrees(groundNormal);
            if (globalInput.keyboard.arrowRight.isDown){
                ax = .7f; // Accelerate right
            }else if (globalInput.keyboard.arrowLeft.isDown){
                ax = -.7f; // Accelerate left
            }else{
                ax = Dot(gameState.rectSpeed, groundDir)*-.1f; // Decelerate.
            }
            gameState.rectSpeed += ax*groundDir;

            float minSpeedY = -30.f;
            bool jumped = false;
            if (grounded){
                if (globalInput.keyboard.arrowUp.isDown){
                    gameState.rectSpeed.y = 10.f; // Jump
                    jumped = true;
                }else{
                    float diff = Abs(AngleDifference(PI/2, AngleOf(groundNormal)));
                    if (diff > PI/4 && gameState.rectSpeed.y > minSpeedY)
                        gameState.rectSpeed.y -= .3f;
                }
            }else if (gameState.rectSpeed.y > minSpeedY){
                gameState.rectSpeed.y -= .4f; // Gravity
            }
            
            // Limit speed
            float maxSpeed = 12.f;
            if (Length(gameState.rectSpeed) > maxSpeed)
                gameState.rectSpeed *= maxSpeed/Length(gameState.rectSpeed);

            float speedLength = Length(gameState.rectSpeed);
            float toMove = speedLength;
            bool bounced = false;
            for(int i = 0; i < 6; i++){ // Artificial iteration limit
                bool collided;
                v2 collisionNormal;
                v2 oldPos = gameState.rectPos;
                v2 DEBUG_test = MoveRectangle(gameState.rectPos, gameState.rectDim/2, V2FromLengthDir(toMove, AngleOf(gameState.rectSpeed)), &collided, &collisionNormal);
                
                gameState.rectPos = MoveRectangle(gameState.rectPos, gameState.rectDim/2, V2FromLengthDir(toMove, AngleOf(gameState.rectSpeed)), &collided, &collisionNormal);
                toMove -= Length(gameState.rectPos - oldPos);
                
                if (!collided)
                    break;

                v2 effectiveSpeed = V2FromLengthDir(speedLength, AngleOf(gameState.rectSpeed));
                if (Dot(effectiveSpeed, collisionNormal) < -9.f && Dot(collisionNormal, V2(0, 1)) < .4f) {
                    // Hit a wall hard: bounce!
                    gameState.rectSpeed = (gameState.rectSpeed - 2*Dot(gameState.rectSpeed, collisionNormal)*collisionNormal)*.4f;
                    bounced = true;
                }else{
                    v2 prevSpeed = gameState.rectSpeed;
                    gameState.rectSpeed = gameState.rectSpeed - Dot(gameState.rectSpeed, collisionNormal)*collisionNormal;

                    if (Abs(collisionNormal.y) > .8f) // The edge is very horizontal
                        if (Abs(Dot(Normalize(prevSpeed), collisionNormal)) > .5f) // and we're going quite perpendicular to the edge
                            break; // Don't slide.

                    if (Abs(Dot(Normalize(prevSpeed), collisionNormal)) > .9f) // Going very perpendicular to the edge
                        break; // Don't slide.
                }
                if (toMove < .1f)
                    break;
            }

            // Don't fly off ascending slopes!

            if (true){ // Set to false to disable sticking to ground.
                // Get the grounded state after moving
                bool newGrounded;
                v2 newGroundNormal;
                MoveRectangle(gameState.rectPos, gameState.rectDim/2, V2(0, -1.f), &newGrounded, &newGroundNormal);
                v2 speedDir = Normalize(gameState.rectSpeed);
                if (!bounced && !jumped && grounded && !newGrounded && Abs(speedDir.y) < .95f){ // If stopped being grounded, and not travelling extremely vertically...
                    // We'll try to move the entity down.
                    float maxDistance = 40.f; // The max distance we'll move down (independent of speed)
                    float maxSlope = 3.f; // This further restricts maxDistance based on the horizontal speed.
                    float deltaX = Abs(gameState.rectPos.x - oldPos.x);
                    float finalMaxDistance = Min(maxDistance, deltaX*maxSlope);

                    bool foundGround;
                    v2 newPos = MoveRectangle(gameState.rectPos, gameState.rectDim/2, V2(0, -finalMaxDistance), &foundGround, &newGroundNormal);
                    if (foundGround && deltaX){
                        float prevY = oldPos.y;
                        bool tooMuchSlope = false;
                        int steps = Ceil(deltaX/2);
                        for(int i = 0; i < steps; i ++){
                            float t = (float)(i + 1)/(float)steps;
                            v2 p = {Lerp(oldPos.x, gameState.rectPos.x, t),
                                    Lerp(oldPos.y, gameState.rectPos.y, t)};
                            finalMaxDistance = Min(maxDistance, t*deltaX*maxSlope);
                            v2 intermediatePos = MoveRectangle(p, gameState.rectDim/2, V2(0, -finalMaxDistance), &foundGround, &newGroundNormal);
                            
                            float dy = intermediatePos.y - prevY;
                            float stepSize = (gameState.rectPos.x - oldPos.x)/(float)steps;
                            float slope = Abs(dy)/stepSize;
                            if (slope > maxSlope){
                                tooMuchSlope = true;
                                break;
                            }
                            prevY = intermediatePos.y;
                        }
                        if (!tooMuchSlope){
                            gameState.rectPos = newPos;
                            gameState.rectSpeed = gameState.rectSpeed - Dot(gameState.rectSpeed, newGroundNormal)*newGroundNormal;
                        }
                    }
                }
            }
        }



        //
        // Render
        //
        glViewport(0,0, windowDim.x, windowDim.y);
    
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glMatrixMode(GL_PROJECTION);
        float a = (windowDim.x ? 2.0f/windowDim.x : 1.f);
        float b = (windowDim.y ? 2.0f/windowDim.y : 1.f);
        float proj[] = {    a,     0,     0,    0,
                            0,     b,     0,    0,
                            0,     0,   1.f,    0,
                         -1.f, -1.0f,     0,  1.f };
        glLoadMatrixf(proj);
    

        glClearColor(.8f, .82f, .7f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        
        //
        // Render world
        //
        
        // Draw walls (fill)
        glBegin(GL_TRIANGLES);
        for(int i = 0; i < gameState.numWalls; i++){
            glColor4f(1.f, .32f, .4f, .8f);
            glVertex2f((gameState.walls[i].p[0].x - viewPos.x)*gameState.cameraScale, (gameState.walls[i].p[0].y - viewPos.y)*gameState.cameraScale);
            glVertex2f((gameState.walls[i].p[1].x - viewPos.x)*gameState.cameraScale, (gameState.walls[i].p[1].y - viewPos.y)*gameState.cameraScale);
            glVertex2f((gameState.walls[i].p[2].x - viewPos.x)*gameState.cameraScale, (gameState.walls[i].p[2].y - viewPos.y)*gameState.cameraScale);
        }       
        glEnd();

        // Draw walls (outline)
        glBegin(GL_LINES);
        for(int i = 0; i < gameState.numWalls; i++){
            wall *w = &gameState.walls[i];
            glColor4f(0.02f, 0.f, 0.1f, .8f);
            glVertex2f((w->p[0].x - viewPos.x)*gameState.cameraScale, (w->p[0].y - viewPos.y)*gameState.cameraScale);
            glVertex2f((w->p[1].x - viewPos.x)*gameState.cameraScale, (w->p[1].y - viewPos.y)*gameState.cameraScale);
            
            glVertex2f((w->p[1].x - viewPos.x)*gameState.cameraScale, (w->p[1].y - viewPos.y)*gameState.cameraScale);
            glVertex2f((w->p[2].x - viewPos.x)*gameState.cameraScale, (w->p[2].y - viewPos.y)*gameState.cameraScale);

            glVertex2f((w->p[2].x - viewPos.x)*gameState.cameraScale, (w->p[2].y - viewPos.y)*gameState.cameraScale);
            glVertex2f((w->p[0].x - viewPos.x)*gameState.cameraScale, (w->p[0].y - viewPos.y)*gameState.cameraScale);
            
            // Normals
            //glColor4f(1.f, 1.f, 0.1f, .8f);
            //float d = 15.f;
            //glVertex2f(((w->p[0].x + w->p[1].x)/2 - viewPos.x)*gameState.cameraScale, ((w->p[0].y + w->p[1].y)/2 - viewPos.y)*gameState.cameraScale);
            //glVertex2f(((w->p[0].x + w->p[1].x)/2 + d*w->normals[0].x - viewPos.x)*gameState.cameraScale, ((w->p[0].y + w->p[1].y)/2 + d*w->normals[0].y - viewPos.y)*gameState.cameraScale);

            //glVertex2f(((w->p[1].x + w->p[2].x)/2 - viewPos.x)*gameState.cameraScale, ((w->p[1].y + w->p[2].y)/2 - viewPos.y)*gameState.cameraScale);
            //glVertex2f(((w->p[1].x + w->p[2].x)/2 + d*w->normals[1].x - viewPos.x)*gameState.cameraScale, ((w->p[1].y + w->p[2].y)/2 + d*w->normals[1].y - viewPos.y)*gameState.cameraScale);

            //glVertex2f(((w->p[2].x + w->p[0].x)/2 - viewPos.x)*gameState.cameraScale, ((w->p[2].y + w->p[0].y)/2 - viewPos.y)*gameState.cameraScale);
            //glVertex2f(((w->p[2].x + w->p[0].x)/2 + d*w->normals[2].x - viewPos.x)*gameState.cameraScale, ((w->p[2].y + w->p[0].y)/2 + d*w->normals[2].y - viewPos.y)*gameState.cameraScale);
        }       
        glEnd();

        
        // Draw circle
        glColor4f(0.70f, 0.40f, .83f, .85f);
        if (CircleWallCollision(gameState.circlePos, gameState.circleRadius))
            glColor4f(1.f, 0, 0, .8f);
        glBegin(GL_TRIANGLES);
        for(int i = 0; i < 32; i++){
            glVertex2f((gameState.circlePos.x - viewPos.x)*gameState.cameraScale, (gameState.circlePos.y - viewPos.y)*gameState.cameraScale);
            glVertex2f((gameState.circlePos.x + gameState.circleRadius*Cos(i*2*PI/32.f) - viewPos.x)*gameState.cameraScale, (gameState.circlePos.y - gameState.circleRadius*Sin(i*2*PI/32.f) - viewPos.y)*gameState.cameraScale);
            glVertex2f((gameState.circlePos.x + gameState.circleRadius*Cos((i + 1)*2*PI/32.f) - viewPos.x)*gameState.cameraScale, (gameState.circlePos.y - gameState.circleRadius*Sin((i + 1)*2*PI/32.f) - viewPos.y)*gameState.cameraScale);
        }
        glEnd();

        // Draw rect
        glColor4f(0.20f, 0.24f, 0.93f, .8f);
        if (RectangleWallCollision(gameState.rectPos - gameState.rectDim/2, gameState.rectPos - gameState.rectDim/2 + gameState.rectDim)){
            glColor4f(.73f, .1f, .1f, .8f);
        }
        //glColor4f(0.f, 0.f, 0.f, .8f);
        glBegin(GL_QUADS);
        glVertex2f((gameState.rectPos.x - gameState.rectDim.x/2 - viewPos.x)*gameState.cameraScale, (gameState.rectPos.y - gameState.rectDim.y/2 - viewPos.y)*gameState.cameraScale);
        glVertex2f((gameState.rectPos.x + gameState.rectDim.x/2 - viewPos.x)*gameState.cameraScale, (gameState.rectPos.y - gameState.rectDim.y/2 - viewPos.y)*gameState.cameraScale);
        glVertex2f((gameState.rectPos.x + gameState.rectDim.x/2 - viewPos.x)*gameState.cameraScale, (gameState.rectPos.y + gameState.rectDim.y/2 - viewPos.y)*gameState.cameraScale);
        glVertex2f((gameState.rectPos.x - gameState.rectDim.x/2 - viewPos.x)*gameState.cameraScale, (gameState.rectPos.y + gameState.rectDim.y/2 - viewPos.y)*gameState.cameraScale);
        glEnd();

        
        // Draw circle speed
        {
            glColor4f(0, 0, 0, .7f);
            glBegin(GL_LINES);
            glVertex2f((gameState.circlePos.x - viewPos.x)*gameState.cameraScale, (gameState.circlePos.y - viewPos.y)*gameState.cameraScale);
            float f = 3.f;
            glVertex2f((gameState.circlePos.x + f*gameState.circleSpeed.x - viewPos.x)*gameState.cameraScale, (gameState.circlePos.y + f*gameState.circleSpeed.y - viewPos.y)*gameState.cameraScale);
            glEnd();
        }
        // Draw rect speed
        {
            glColor4f(0, 0, 0, .7f);
            glBegin(GL_LINES);
            glVertex2f((gameState.rectPos.x - viewPos.x)*gameState.cameraScale, (gameState.rectPos.y - viewPos.y)*gameState.cameraScale);
            float f = 3.f;
            glVertex2f((gameState.rectPos.x + f*gameState.rectSpeed.x - viewPos.x)*gameState.cameraScale, (gameState.rectPos.y + f*gameState.rectSpeed.y - viewPos.y)*gameState.cameraScale);
            glEnd();
        }

        // Draw placing wall
        if (gameState.numWalls < ArrayCount(gameState.walls)){
            wall *w = &gameState.walls[gameState.numWalls];
            if (placingWallPointIndex > 0 || globalInput.keyboard.control.isDown){
                // Highlight new points
                int numPointsToHighlight = placingWallPointIndex + 1;
                for(int i = 0; i < numPointsToHighlight; i++){
                    glBegin(GL_QUADS);
                    glColor4f(1.f, .7f, .0f, .8f);
                    float r = 10.f;
                    glVertex2f((w->p[i].x - viewPos.x)*gameState.cameraScale - r, (w->p[i].y - viewPos.y)*gameState.cameraScale - r);
                    glVertex2f((w->p[i].x - viewPos.x)*gameState.cameraScale + r, (w->p[i].y - viewPos.y)*gameState.cameraScale - r);
                    glVertex2f((w->p[i].x - viewPos.x)*gameState.cameraScale + r, (w->p[i].y - viewPos.y)*gameState.cameraScale + r);
                    glVertex2f((w->p[i].x - viewPos.x)*gameState.cameraScale - r, (w->p[i].y - viewPos.y)*gameState.cameraScale + r);
                    glEnd();
                }
                // Show line between the first two points.
                if (placingWallPointIndex == 1){
                    glBegin(GL_LINES);
                    glColor4f(.9f, .7f, .1f, .8f);
                    glVertex2f((w->p[0].x - viewPos.x)*gameState.cameraScale, (w->p[0].y - viewPos.y)*gameState.cameraScale);
                    glVertex2f((w->p[1].x - viewPos.x)*gameState.cameraScale, (w->p[1].y - viewPos.y)*gameState.cameraScale);
                    glEnd();
                }
                // Show triangle when placing the third point.
                if (placingWallPointIndex == 2){
                    glBegin(GL_TRIANGLES);
                    glColor4f(1.f, .7f, .2f, .5f);
                    glVertex2f((w->p[0].x - viewPos.x)*gameState.cameraScale, (w->p[0].y - viewPos.y)*gameState.cameraScale);
                    glVertex2f((w->p[1].x - viewPos.x)*gameState.cameraScale, (w->p[1].y - viewPos.y)*gameState.cameraScale);
                    glVertex2f((w->p[2].x - viewPos.x)*gameState.cameraScale, (w->p[2].y - viewPos.y)*gameState.cameraScale);
                    glEnd();
                }
            }
        }


        glFlush();
        SwapBuffers(dc);



        //
        // Sleep to render at 60 FPS
        //
        while(true){
            LARGE_INTEGER newFrameTime = GetCurrentTimeCounter();
            float timeElapsed = GetSecondsElapsed(lastFrameTime, newFrameTime);
            if (timeElapsed > 1/60.f){
                lastFrameTime = newFrameTime;
                break;
            }
            if (1/60.f - timeElapsed > 0.005f){
                Sleep(1);
            }
        }

        // Reset button input.
        for(int i = 0; i < ArrayCount(globalInput.keyboard.asArray); i++){
            globalInput.keyboard.asArray[i].transitionCount = 0;
        }
        for(int i = 0; i < ArrayCount(globalInput.mouseButtons); i++){
            globalInput.mouseButtons[i].transitionCount = 0;
        }
    }

    if (CREATE_CONSOLE){
        FreeConsole();
    }

    wglMakeCurrent(NULL, NULL);
    ReleaseDC(window, dc);
    wglDeleteContext(rc);
    DestroyWindow(window);

    return 0;
}
