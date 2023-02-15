# How to program a simple vectorial 2D platformer

In this tutorial, I will show you how to program a vectorial 2D platform game from the ground up. The code is in straightforward, procedural C++, but all the concepts are language-agnostic. The tutorial focuses on the fundamentals of collision detection and the details that make this simple platformer engine work.

The prerequisites for this tutorial are C/C++ knowledge and high school level math. The public project only targets Windows, but I made a Wasm version that I [put up on a website](https://pered.itch.io/platformer-tutorial-demo). 

TLDR: [See the central function in the source code](https://github.com/Pere001/2d-platformer-tutorial-2023/blob/main/src/main.cpp#L717). [Play a web version of the final demo](https://pered.itch.io/platformer-tutorial-demo).

__Table of Contents__

* [Intro](#intro)
  * [The Problem](#the-problem)
  * [Overview of my Method](#overview-of-my-method)
* [Part 0 - Basic Math](#part-0---basic-math)
  * [Trigonometric Functions](#trigonometric-functions)
  * [Dot and Cross Products](#dot-and-cross-products)
* [Part 1 - Collision Functions](#part-1---collision-functions)
  * [Computing Normals](#computing-normals)
  * [Rectangle vs. Triangle](#rectangle-vs-triangle)
  * [Circle vs. Triangle](#circle-vs-triangle)
* [Part 2 - Move in Incrementally Smaller Steps](#part-2---move-in-incrementally-smaller-steps)
* [Part 3 - Collision Normals](#part-3---collision-normals)
  * [Colision Normal: Circle vs. Triangle](#collision-normal-circle-vs-triangle)
  * [Collision Normal: Rectangle vs. Triangle](#collision-normal-rectangle-vs-triangle)
* [Part 4 - Using the Normal](#part-4---using-the-normal)
  * [Sticking to the Ground](#sticking-to-the-ground)
* [Part 5 - Optimization](#part-5---optimization)
* [Closing Words](#closing-words)
  * [Further Reading](#further-reading)


# Intro

![Demo](/img/demo.gif)

## The Problem

Every platformer game has to handle movement and collisions. Therefore, we must design a system that can query the physical situation of an entity (e.g. grounded state, raycasts, etc.) and apply movement in a way that conforms predictably to the terrain.

This can be a daunting task because you need to solve many open problems in a highly coordinated manner. 

For example, the entity state might be useful in the process of resolving the collisions, which might in turn change the entity state and inform the gameplay code, whose needs will determine the kind of information we need to query about the surroundings of an entity, which will then determine what entity and terrain state we need to store, etc. Each part imposes constraints and demmands upon the others, so it's hard to know where to start.

I can imagine three kinds of approaches to designing the core of a platformer:

1. Depend on a physics engine. A physics engine is a complex system that follows the equations of Newtonian mechanics to systematically apply realistic movement to objects that might have different properties and constraints. It's what you need for games like _Angry Birds_. If your platformer features objects with this kind of realistic physics, you could make your player character a physical object as well, and control it indirectly through the interface of the physics engine (e.g. applying forces). This is probably what they do in games like _Limbo_ and _Little Big Planet_. 
    
    In this tutorial, we won't cover physics engines. We will focus on the nonrealistic types of platformers. These can still emulate some effects from real life physics (for example, all platformers have some interpretation of gravity), but they're not structured to solve complex physical interactions of many solid objects.

2. What I would call "the railroad method". If your entity is grounded, keep track of the point it's grounded on, associated with an edge of a specific wall. This requires the different edges of a continuous floor to be linked by some metadata. In this method, it's OK if the entity slightly intersects some walls while it's grounded, because its walking movement is bound to the chain of edges like on a railroad.

    I think this is the kind of method they use in games such as _Rayman Origins_ and _Braid_. This is probably more complex than the method we'll use and also more specialized, since it requires some coupling with the player movement mechanics and some specific formatting of the terrain data.

3. My solution to the problem is a guarantee of separation. We ensure no entity ever intersects any wall. We check for collision before attempting to move an entity to a new place. If the target position collides, we make smaller and smaller safe steps on the segment of the speed until the steps are so small that it looks like the entity is touching a wall. However, it is always slightly separated. 

    To avoid all these iterations, we could find the tangent position directly, add a little margin to it, and move the entity there. But, for this tutorial, we'll just use the stupid iterative approach to keep things simple. With this method, there's no "terrain metadata" at all. We simply use the raw geometry and the assumption of separation to find normals that we'll use to guide the entity's movement.

I'm sure there are other possible options, and games could also use a mixed approach.

## Overview of my Method

Here's a summary of the method (don't worry if you don't totally get it right now). We try to move an entity by its speed. If it intersects any wall, we instead try to move it by a fraction of the speed. Through "convergent" trial and error tests, we find the highest fraction of the speed that we can add to the entity's position without intersecting any wall, and we move the entity there.

Then, we use some math to find the first wall to collide and its collision normal. Finally, we use the normal to update the speed, with the possibility of applying bounce or slide.

Additionally, there's a trick we do to avoid being stuck because of precision issues when moving parallel and very close to a wall's edge. When we're testing positions we will also test if "shifting" the entity by a subpixel amount perpendicular to the speed allows it to go farther, and if so move to that farthest, shifted position.

This method has a downside. It's rather costly because it does many iterations of collision checks against all candidate walls. On top of that, the trick that avoids getting stuck when moving tangent to edges potentially triples the number of iterations.

That said, this method is simple and modular. You can easily improve or simplify different aspects to fit your needs. You can add support for a different shape of entity or wall by supplying the respective intersection function and the code for the collision normal. This approach is also robust to very rough terrain data. Your walls can be paper-thin. They can even intersect one another. You will still get the desired behavior.

With the math described in this tutorial, you can improve the movement code so it doesn't do the repeated iterations of intersection checks to solve collisions. Instead, you would directly find the farthest position that doesn't collide by calculating the time of impact and advancing only by that time. (Honestly, I'd do that if I were remaking this tutorial from scratch, now that I know better. But it's also kind of interesting to experiment with these more wacky ideas, and it does keep the tutorial simpler.)

The structure of this tutorial follows the steps I'd recommend as the implementation arc:

0. [Basic Math](#part-0---basic-math): Of course, you will need your math functions, but most important is an excellent intuition of the dot and cross products.
1. [Collision Functions](#part-1---collision-functions): We will start by making the collision functions. After this, you can test them by placing a bunch of walls and dragging a shape with the mouse that changes color when it intersects a wall.
2. [Move in Incrementally Smaller Steps](#part-2---move-in-incrementally-smaller-steps): Next, we will implement the iterative testing I described, to move an entity up to the edge of a wall without intersecting it. You can test this by creating an entity that is movable by keyboard input.
3. [Collision Normals](#part-3---collision-normals): We will calculate the collision normals. You can test this by drawing a line in the direction of the last normal.
4. [Using the Normal](#part-4---using-the-normal): We will use the normal to update the player's speed, and we will finally have a moving player.
5. [Optimization](#part-5---optimization): We will finish by talking about optimizations and improvements you can make to this project.

For entities, I chose to implement a circle and a rectangle. For walls, I chose triangles. With this approach, it is easy to add more entity shapes, and it is trivial to change the walls from triangles to general convex polygons. You just have to change the hard-coded `3` in a bunch of `for` loops.

Without further ado, let's get on with the tutorial.


# Part 0 - Basic Math

Before we start, let's clarify the mathematical functions and conventions we will use. Directions will be implicitly in counterclockwise radians starting from the right. We will define a 2D vector struct (`v2`) and overload all the operators corresponding to its mathematical operations (I assume you understand basic vector math).

Our coordinate system will be positive y up. We will be using `float`s for our positions, but note that you could also use `int`s with very little change to the code, as I do in my own game.

```c++
struct v2{
    float x;
    float y;
};
```

The rest of this section is mainly intended for beginners.

## Trigonometric Functions

Here are the basic trigonometric functions we will use (in pseudocode).

- `Length(v) = SquareRoot(v.x*v.x + v.y*v.y)`
- `LengthSqr(v) = v.x*v.x + v.y*v.y` It’s common to use the length squared as an optimization. When comparing distances with <, >, <= or >= the result is the same as comparing those distances squared, so sometimes we can save ourselves the square root.
- `AngleOf(v) = Atan2(v.y, v.x)` (Since atan2 is undefined for zero vectors, we will also check if v.x and v.y are 0 and return 0 in that case.)
- `V2FromLengthDir(length, dir) = {Cos(dir)*length, Sin(dir)*length}` Constructs a `v2` from a length and an angle.
- `Normalize(v) = V2FromLengthDir(1, AngleOf(v))` I like this option better than the common alternative `v/Length(v)` because it allows the input to be a zero vector (since it avoids the division by zero). This will reduce the amount of edge cases needed in our code.
- `Rotate90Degrees(v) = {-v.y, v.x}`
- `RotateMinus90Degrees(v) = {v.y, -v.x}`

## Dot and Cross Products

- `Dot(a, b) = a.x*b.x + a.y*b.y`
- `Cross(a, b) = a.x*b.y - a.y*b.x`

The dot product between 2 vectors, if one of them is unit-length (normalized), is a projection of the non-unit-length vector into the direction of the unit-length one. That is, what (signed) distance along the axis of the unit-length vector the other vector travels.

![Dot product](/img/0_0_dot_product.jpg)

The cross product helps us compare the angle between two vectors. If the two vectors are parallel, the cross product will be 0. If both vectors are normalized, and they are perpendicular to each other, the cross product will be 1 (if the angle from a to b is 90°) or -1 (if the angle is 270°).

Most importantly, if the second vector is less than pi radians (180°) counterclockwise from the first, the cross product will be positive, and if it’s over pi radians it wil be negative.

![Cross product](/img/0_1_cross_product.jpg)

[Check out this Desmos graphic to see this interactively](https://www.desmos.com/calculator/foxjolfgkq). I admit I often come to play with it when I need to refresh my understanding of the two products. It is crucial that you have a basic intuition about the dot and cross products if you are to follow the rest of this tutorial.

I must note that in mathematics the "cross product" is technically only defined for 3D vectors. In 3D it’s a vector perpendicular to vectors a and b, and its formula is: `a × b = (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)`. 

As you can see, the Z component matches the definition of cross product that we use. That’s because in 2D it’s common to use the term “cross product”, as we do, to refer to the Z component of the actual cross product of your 2D vectors extended to 3D. In some code bases they also call this "the determinant”, which is an operation you can do on a matrix that also happens to match our 2D cross product in the case of 2x2 matrices. So yeah, we'll call it the cross product as many people do, but be aware that might not be perfect math language.


# Part 1 - Collision Functions 

## Computing Normals

We always store the vertices of our walls in the same directional order, in our case counterclockwise (CCW), with the normals pointing outwards. We will precompute the normals and store them with the wall data to avoid computing them every time we need them (which will be often).

```c++
struct wall{
    v2 p[3]; // Vertices.
    v2 normals[3]; // normals[0] is the normal of the edge from p[0] to p[1], and so forth.
};
```

To calculate the normal of an edge we just do:

```c++
v2 edgeDir = Normalize(p1 - p0);
v2 n = {edgeDir.y, -edgeDir.x}; // Rotates edgeDir -90 degrees CCW
```

To ensure three points are CCW we can simply check the cross product of the first and second edges:

```c++
bool PointsAreCW(v2 p0, v2 p1, v2 p2){
    return (Cross(p1 - p0, p2 - p1) < 0);
}
```

Just for completion, if you had a general polygon instead of a triangle, you could use this slightly more complex method. Sum the angle of the direction of each edge relative to the direction of the previous edge. If the points make a full circle CCW, the sum will be `2*PI` and if they make a full circle CW the sum will be `-2*PI`.

```c++
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
```

You could flip the order of the points of a wall if this returned true. Or you could just load the walls from valid data.

## Rectangle vs Triangle

To check whether a rectangle and a triangle intersect, we will use the Separating Axis Theorem (SAT), which is probably the most common method for finding collisions between convex polygons. If we can find one line that separates all of the rectangle's points from all of the triangle's, we will know they do not intersect.

The SAT says that if two convex polygons are separated, one of the lines that separate them will be in the direction of one of the edges of our polygons. This means we only have to check as many directions as edges there are in both polygons.

To test for separation in the direction of a given edge, we can project all the points of the other shape, relative to a point on the edge, onto the normal of that edge. If all the projections are greater than zero, it means a line in the direction of the edge can separate the walls. We will do this for every edge. If at least one edge finds separation, the shapes do not intersect. If no edges find separation, the shapes do intersect.

I don’t count tangential shapes as colliding, but you could. In our case I don’t think it matters.

```c++
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
```

If you got lost, I suggest thinking by making some drawings, or just check out my diagram:

![Rect vs tri](/img/1_1_rectangle_vs_triangle_1.jpg)

This illustrates the check of an edge that finds separation. We consider all points relative to a vertex of the edge we are checking. All the projections of the triangle vertices onto the normal of that edge will be <= 0, because the normal points outwards, and no point on the triangle can be in front of the normal. If the minimum of the projections of the rectangle is >= 0, it means there is separation.

## Circle vs Triangle

This one is fun. First, we will project the center of the circle onto the normals of the triangle’s edges. This is similar to what we did with the rectangle, and it gives us the signed distance of the center to the line of each of the edges. If that distance is greater than r in all of them, there is no collision.

But this is not enough. Were we to rely solely on this test, we would see false positives near the triangle's vertices. In essence, we would be checking a collision between the circle's center and a larger triangle with each edge at distance r from our original triangle. This is not equivalent to the collision between our circle and the original triangle.

![Circle vs tri](/img/1_2_circle_vs_triangle_1.jpg)

A collision of the center point against the larger triangle with rounded corners IS equivalent.

![Circle vs tri](/img/1_2_circle_vs_triangle_2.jpg)

So we’ll break the collision down into 3 shapes:

![Circle vs tri](/img/1_2_circle_vs_triangle_3.jpg)

1. If the center is inside the (original) triangle, that’s a collision.
2. If the center is inside any of the rectangles created by extruding the triangle’s edges by r (circle’s radius), that’s a collision.
3. If the center is inside a circle of radius r with center on any of the triangle vertices, that’s a collision.

Otherwise there’s no collision.

```c++
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
```
The first loop does this, with respect to the position of the center:

![Circle vs tri](/img/1_2_circle_vs_triangle_4.jpg)

The second loop just returns true if the center is inside those vertex circles (which is the same as checking if the vertices are inside the original circle).


# Part 2 - Move in Incrementally Smaller Steps

Now that we have the collision functions, let’s apply some velocity. We attempt to move the entity to the new position (old position + speed), and if we don’t collide, we’re done. If there is a collision at the new position, we'll check multiple points on the segment between the old and new positions. We will then keep the farthest point that doesn't collide with any wall. If all of them collide, we will keep the old position.

How do we decide which intermediate points to check? We could move in equally small steps toward the new position until we collide, but there is a more efficient way. We will cut the step size in half at every iteration. This way, for the same level of precision, we will only have to do O(log(n)) checks instead of O(n) (see [big O notation](https://en.wikipedia.org/wiki/Big_O_notation)).

![Incrementally smaller steps](/img/2_1_incrementally_smaller_steps.jpg)

First, the step size is the length of speed. We advance toward pos + speed by the step size. If we collide, we half the step size and move backwards. If that collides, we half the step and move backwards again. If we don’t collide, we half the step and move forwards. Etc. We’ll stop when we reach an arbitrary limit on the step size or on the number of iterations. 0.05 pixels for the minimum step length and 15 maximum steps seems reasonable enough. You can tune these for different kinds of entities to save some performance.

Here is the first version of the `MoveCircle()` function:

```c++
v2 MoveCircle(v2 pos, float r, v2 speed){
    v2 newPos = pos;
    float stepSize = 1.0f;
    float speedFactor = 1.0f;
    for(int steps = 0; steps < 15; steps++){
        v2 p = pos + speedFactor*speed;
        stepSize *= 0.5f;
        if (CircleWallCollision(p, r)){ // (This calls CircleWallCollision on all existing walls)
            speedFactor -= stepSize;
        }else{
            speedFactor += stepSize;
            newPos = p;
            if (steps == 0)
                break; // No collisions found on the first step.
        }
        if (stepSize*Length(speed) <= .05f)
            break; // Limit step size.
    }
    return newPos;
} 
```

I kept this code simple because we will expand upon it later.
- After moving the circle we will want to compute the collision normals if it collided.
- We will add small perpendicular shifts to facilitate moving tangent to edges. This solves the problem of entities getting stuck when sliding along edges. Because of precision issues, if you move a shape along an edge, it will find collision at some positions and not others.

    Our solution shifts the entity by a subpixel amount perpendicular to its speed when a collision is found. If the shifted position doesn't collide, we will move there and count that as a non-collision (so, the next step will be forward rather than backward). You can see the implementation of this in [the final version of this function](https://github.com/Pere001/2d-platformer-tutorial-2023/blob/main/src/main.cpp#L574).
- There will be an equivalent function for the rectangle, `MoveRectangle()`.
- To avoid tunnelling (fast entities going through thin walls) we could have a maximum step length by adding an outer loop.


# Part 3 - Collision Normals

If we can't move the entity by the whole speed, we will want to find the collision normal of the first wall the entity collides with. Sometimes this is the normal of one of the wall’s edges, and sometimes it will be something else (imagine the entity colliding against the protruding edge of a wall, like a basket ball hitting the corner of a table).

We need the collision normal because it helps us calculate the bouncing or sliding reaction we will see in the next section, along with other possible mechanics.

As a recap, these are the data we'll use to calculate the collision normal after having moved the entity and detected a collision:
- Entity speed.
- New entity position, or the position of the last non-colliding step. It shouldn't collide with any wall.
- Last colliding position, or the position of the entity at the last colliding step. This should be farther than the new position along the direction of the speed, and infinitesimally close to it.
- The walls. We will use the ones that intersect the entity in its last colliding position.

## Collision Normal: Circle vs Triangle

The following is a simple method to get the circle's collision normal. Find the closest point on the closest edge of the closest wall (to our new position). The collision normal will be the direction from that point to the center of the circle.

![Circle vs triangle normal](/img/3_1_circle_vs_triangle_normal_1.jpg)

The blue arrow is the circle’s speed. The black segment is the closest edge. The black dot is the closest point on that edge. The black arrow is the resulting collision normal.

This distance-based method works very well because the circle has previously been moved very close to the wall it will hit. However, it's still possible that the edge that is closest is not the one it would hit first. We will solve this by taking the velocity's direction into account.

Instead of comparing the edges by their distance to the circle, we will compare them by the time it takes the circle to hit them if it keeps moving by its speed.
To compute this time, we will simply divide the distance by the dot product of the velocity and the direction from the circle to the point. This will prioritize edges that the circle is moving towards. The units also make sense, because if you divide distance by speed, you get time.

But it's probably better understood geometrically:

![Circle vs triangle normal](/img/3_1_circle_vs_triangle_normal_2.jpg)

Now, if we take the lowest time, we correctly find the wall that collides first. The lines representing "time" are actually showing time\*speed, or the distance travelled until collision.

Here's a more detailed proof: 

![Circle vs triangle normal](/img/3_1_circle_vs_triangle_normal_3.jpg)

As you can see, the dot product `Dot(speed, -n)` projects the speed onto `-n` to give us the amount that `speed` is moving towards our edge. That is, the distance the circle will move along the direction of the negative normal each second. Or the distance our circle will move towards our edge each second.

We want to know how much time it will take to travel the distance `d` along the direction of the normal. Since both are along the same direction, we can divide the distance we want to travel by the projected speed, and we will get the time it takes. It's not that complicated when you break it down.

Note that `n` is not really taken from the edge's normal. To find it, you compute the direction from the closest point on the edge, to the center of the circle. This allows us to handle the cases where the circle collides a pointy vertex, using exactly the same code as when it collides the surface of an edge.

Here’s the code. Imagine it’s inside the previous `MoveCircle()` function, after the loop. `walls` is an array of all existing walls, sized by `numWalls`. Imagine we then return the `bestNormal` along with the new position and whether it collided. You can also see [the whole function on GitHub](https://github.com/Pere001/2d-platformer-tutorial-2023/blob/main/src/main.cpp#L574).
```c++
v2 bestNormal = -Normalize(speed); // Default value in case we don't find any.
float bestTime = MAX_FLOAT; // Lowest time until collision, to find the edge that would collide first.
for(int wallIndex = 0; wallIndex < numWalls; wallIndex++){
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
        if (time < bestTime && Dot(n, speed) <= 0){ // Ignore normals that face away from speed (technically possible with shifts)
            bestTime = time;
            bestNormal = n;
        }
    }
}
```

## Collision Normal: Rectangle vs Triangle

For the rectangle's collision normal, we will borrow some techniques we used with the circle. We will iterate each wall that collides the rectangle in its last colliding position to find the wall that collides first.

For each wall, we will iterate every edge of that wall and every edge of the rectangle, treating them the same. For every edge, we will project onto its normal all the points of the other shape. Of these projected points, we'll take the lowest. This will be the "projected distance" of this edge in respect to the other shape. From this "projected distance", we will calculate the "time until impact" of the edge by dividing by the dot product between the speed and the normal (the same way we did with the circle). As I have said, we will do this with all edges of both shapes. The edge with the highest "time until impact" will provide the "time until impact" of a given wall and its collision normal. Why the highest? The best way to understand this is with some visual examples:

![Circle vs triangle normal](/img/3_2_rectangle_vs_triangle_normal_1.jpg)

The grey lines show the "projected distance" of edges B and C. The black lines represent the "time until collision" (a, b, and c) of edges A, B, and C. It is clear to the eye that the correct time until impact is a, which is the highest. It's also clear that the collision normal will be the normal of edge A.

So, taking the highest "time until impact", we will find the "time until impact" for a wall and its collision normal. We will do this for all candidate walls and the winner will be the wall with the lowest time. The collision normal of this wall will be the final collision normal. If the collision normal was taken from an edge of the rectangle, we will negate it, because we want it from the point of view of the rectangle.

If we stop here, this will give us the expected normal in most situations. However, it might cause problems with geometry such as the following:

![Circle vs triangle normal](/img/3_2_rectangle_vs_triangle_normal_3.jpg)

In this example, edge B is clearly in front of edge C, covering it fully, so we would expect the normal of B to win over the normal of C. The collision normal of A is the same as B's, so their respective times, a and b, are equal. Therefor, the upper left wall will select either A or B as the winner edge, but either way the collision normal will be the left direction.

The lower right wall, however, will pick either A or C as the winner edge, since c also happens to be the same as a and b. Its collision normal can be either the left direction or the normal of C. Between the two walls, any of them can win, because we know their "time until impact" will be the same. So, if the lower right wall wins, and it had selected C as its colliding edge, the normal of C will mistakenly be chosen as the final collision normal. This could cause unexpected behaviour, like our player bouncing upwards when it runs into this wall.

To solve this, we must understand which situations can cause this problem to happen. It wouldn't happen if the rectangle was slightly higher up, because the line of C would be farther, c would be greater, and we would correctly collide with edge B. It also wouldn't happen if the rectangle was a bit lower down, because the line of C would be closer, c would be lower than a, and the lower right wall would select A as the colliding edge.

It can only happen if we hit such a point with perfect alignment. Most of the time this happens when there is a floor at the level of the point that leads the rectangle's lower corner into it. So, our solution will need to tackle this case.

When two walls have the same "time until impact", we will break the tie by choosing the wall with a negated normal more aligned with the velocity of the rectangle.

- When the rectangle is walking horizontally to the right, like the example above, the edge with the normal that faces leftmost will win the tie. This will choose the correct edge, B, and it also handles the cases where you're moving uphill.

- In the cases where the rectangle is moving downhill, C wins this preference. However, the rectangle won't get to hit that point of conflict if the floor is connected to it, because it will be stopped by its own hitbox hitting B on a higher position.

Long story short, this simple way to break the tie (just a dot product) solves the problem in the disproportionately most likely situation for it to occur. The unaddressed cases very improbable.

For the glitch to arise now, a rectangle would have to jump into a wall that contains such a conflictive vertex in the perfect trajectory that properly aligned a corner of the rectangle with the conflictive vertex, and then the wrong edge would have to win its 1/4 chance of being selected. You could do something more complex to address these implausible cases. In practice, this will be more than enough for most games.

The equality check for whether there's a tie will be done with an epsilon (i.e. check that the two values are __almost__ the same), since the two values will be computed by different routes, so even if geometrically they should be the same, they might be slightly different.

That's all we need to do to find the correct normal. Here's the code (this goes at the end of [`MoveRectangle()`](https://github.com/Pere001/2d-platformer-tutorial-2023/blob/main/src/main.cpp#L717)):

```c++
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

```


# Part 4 - Using the Normal

Now that the core of the engine is done, we can use the collision normal to update the movement of the player and other entities. There are many ways we could use this information to customize the movement: limits on walkable steepness, sticking to the ground or not, bouncing or not, the way speed is accumulated when going down curved ramps, change of acceleration depending on slope, etc.

We also have all the degrees of freedom of platformers in general: acceleration, speed, maximum velocity, types of jump archs, jump cancelling, "coyote jump", stepping over small obstacles, friction, etc.

I will show how to implement smooth, sliding movement; sticking to the ground; and bouncing. I'll keep the code simple, but note that you could add a bunch of quirks to affect the feel of the movement. I'll be showing the update code for the rectangle, but I've done a similar thing for the circle as well.

First, we get the grounded state and ground normal by calling `MoveRectangle()` as if we were trying to move 1 pixel downwards.

For reference, here's the prototype of that function: `v2 MoveRectangle(v2 pos, v2 halfDim, v2 speed, bool *outCollided, v2 *outCollisionNormal)`.

We discard the resulting position so the rectangle doesn't actually move yet. We will use the ground normal later on. It helps us determine the direction of the acceleration caused by horizontal input. For this reason, if we're not grounded, we will set the ground normal to point up, so that when we are in the air the horizontal input just causes purely horizontal acceleration. The grounded state will be used to determine if the player can jump or not.

```c++
bool grounded;
v2 groundNormal;
// We just call this to get the grounded state
MoveRectangle(gameState.rectPos, gameState.rectDim/2, V2(0, -1.f), &grounded, &groundNormal);
if (!grounded){
    groundNormal = V2(0, 1.f);
}
```

Next, we query the input to accelerate the entity "horizontally". It's not really horizontal because we apply the acceleration in the direction of the ground so that the entity will "stick" to the ground.

```c++
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
```

Then, we apply vertical acceleration. We will only apply gravity if the player is not grounded. If the player is grounded, we will let it jump. If it doesn't jump and it is grounded on a steep enough incline, we will also add some gravity, so it will naturally slide down the slope. This steepness threshold is arbritrary, and you can tweak it to fit your game (as all the other magic numbers you see around here). Finally, we limit the speed length so that the player doesn't accelerate too much.

```c++
float minSpeedY = -30.f;
if (grounded){
    if (globalInput.keyboard.arrowUp.isDown){
        gameState.rectSpeed.y = 10.f; // Jump
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
```

Now that we have updated the speed, we will call `MoveRectangle()` to update the position of the entity. If there was a collision, we will use the collision normal to reflect the speed, causing a bounce. We can also project the speed onto the collision edge, causing a slide. The following diagram shows how slide or bounce is achieved by updating the speed (s) based on the collision normal (n).

![Slide and bounce](/img/4_0_slide_bounce.jpg)

The decision of whether to slide or bounce is made from a handful of arbitrary heuristics tuned to my taste: namely, the projection of the speed onto the collision normal, and the direction of the collision normal (i.e. steepness).

This is enough to make the entity react properly to the terrain, but it will make the entity's movement stuttery every time there is a collision, because in these cases the entity moves less than the speed. The following diagram depicts the position of an entity in consecutive frames as it walks on some irregular ground:

![Smoothness](/img/4_0_smoothness_1.jpg)

You can see that where the slope begins, the entity doesn't advance by the whole speed. To solve this, we will do multiple iterations. We will keep track of how much of the original speed we have advanced. If we slid and didn't advance enough, we will do another iteration, now advancing in the new direction, until we've consumed all the speed, or an iteration limit is reached:

![Smoothness](/img/4_0_smoothness_2.jpg)

So that's what this `for` loop does:

```c++
float speedLength = Length(gameState.rectSpeed);
float toMove = speedLength;
for(int i = 0; i < 6; i++){ // Artificial iteration limit
    bool collided;
    v2 collisionNormal;
    v2 oldPos = gameState.rectPos;
    gameState.rectPos = MoveRectangle(gameState.rectPos, gameState.rectDim/2, V2FromLengthDir(toMove, AngleOf(gameState.rectSpeed)), &collided, &collisionNormal);
    toMove -= Length(gameState.rectPos - oldPos);
    
    if (!collided)
        break;

    v2 effectiveSpeed = V2FromLengthDir(speedLength, AngleOf(gameState.rectSpeed));
    if (Dot(effectiveSpeed, collisionNormal) < -9.f && Dot(collisionNormal, V2(0, 1)) < .4f) {
        // Hit a wall hard: bounce!
        gameState.rectSpeed = (gameState.rectSpeed - 2*Dot(gameState.rectSpeed, collisionNormal)*collisionNormal)*.4f;
    }else{
        v2 prevSpeed = gameState.rectSpeed;
        gameState.rectSpeed = gameState.rectSpeed - Dot(gameState.rectSpeed, collisionNormal)*collisionNormal;

        if (Abs(Dot(collisionNormal, V2(0, 1))) > .8f) // The edge is very horizontal
            if (Abs(Dot(Normalize(prevSpeed), collisionNormal)) > .5f) // and we're going quite perpendicular to the edge
                break; // Don't slide.

        if (Abs(Dot(Normalize(prevSpeed), collisionNormal)) > .9f) // Going very perpendicular to the edge
            break; // Don't slide.
    }
    if (toMove < .1f)
        break;
}
```

And there we have it, a working platformer:

![No sticking](/img/no_sticking.gif)

As you can see, when we're on the slope, the rectangle moves in the direction of the floor. However, the rectangle doesn't stick to the ground when there's a decreasing change in slope, and we fly off. This might be the behaviour you want, but just for completion, I'll show a quick way to implement truly sticking to the ground.

## Sticking to the Ground

After moving the entity, we will check its grounded state again. If it was grounded before moving, and it's not grounded after moving, it has flied off. That is, assuming it didn't bounce and that it didn't jump. So, we will add the simple variables `jumped` and `bounced` that we will set to true in the parts of the previous code that handle jumping and bouncing respectively. The entity could also have stopped being grounded because it slid into a vertical edge. We don't want to pull it down in that case, so we will also check the direction of the speed, and if it's highly vertical we won't count that as having flied off.

Now that we have a way to detect whether the entity has flied off, we will try to pull it down if it has. We will define a maximum change in slope (i.e. steepness) we can correct, in our case 3. Then we will find the maximum distance we will allow pulling it down, based on that slope and the distance it travelled horizontally. We will use `MoveRectangle()` yet again to attempt to move the entity downward by this max distance. If we detect a collision, it will mean the entity has moved down into a grounded position. We can take this as the new position. If there was no collision, it means no ground was found, so we can leave the entity as it was, letting it fly off. 

Depending on the terrain this might be enough. But in situations such as the following, it will snap to the ground when it shouldn't:

![Sticking to the ground 1](/img/4_1_sticking_to_the_ground_1.jpg)

In this case, the entity shouldn't stick to the ground because there is an abrupt change in the floor line. We will detect these abrupt changes by iterating from the old position (1 in the diagram) to the new (2) by small steps, and each step finding the distance to the ground by calling `MoveRectangle()` as we did before. If the difference between the ground distance in two consecutive steps is greater than the maximum slope we allow times the step size, we will not pull the entity down. For this loop, we will need the old position (1), so before the movement code we will store the position of the entity into a variable `oldPos`.

This is the code responsible for sticking to the ground, after the movement code:

```c++
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
```

You can see how the entity now sticks to the ground:

![Sticking](/img/sticking.gif)

Note that this is an inefficient manner of sticking to the ground, because of the amount of iteration. But if you don't have a ton of entities running this and constantly flying off the ground, it will work fine.

Also note that this can cause a small irregularity in the speed, because by shifting the entity down, the distance travelled in that frame changes. However, this is probably imperceptible enough that you don't need to worry about it.


# Part 5 - Optimization

We won't be doing any optimization in this article, but I will explain what optimizations you could do if you wanted to make a real game out of this.

First of all, you want to avoid iterating all existing walls when you check for collision. Instead of checking all walls, `MoveCircle()`/`MoveRectanlge()` would take in a list of walls to check against. You could have all your terrain organized into a grid of rectangular chunks, big enough so that no entity can step on walls from 3 different chunks in the same axis in the same frame.

Each wall could be assigned to a chunk. Before calling the `MoveCircle()`/`MoveRectangle()` functions, you make a rectangular bounding box fully containing the entity in all the potential final positions after being moved (you can take the rectangle defined by the speed vector and expand it to account for the entity's dimension, the potential shift, and a tiny margin for precision issues). 

Next, you iterate all walls in the 4 chunks closest to the bounding box and cache those walls that intersect it. After that, you can call `MoveCircle()`/`MoveRectangle()` passing the list of cached walls. If you have a loop that calls the Move functions multiple times for the same entity, as we do in our player movement code, an additional broader caching step might further reduce the number of collision checks. 

Here are some other improvements you might want to make to `MoveCircle()`/`MoveRectangle()`:
- Allow passing the minimum step size and maximum number of iterations as parameters, so you can configure different entities with different precision levels. Also, the ground checks don't need intermediate iterations, although ideally that would have its own code instead of reusing the Move functions.
- Allow passing parameters that can configure the shift amount and let you disable the shifts, again, so you can configure entities differently.
- You could return whether it shifted, if you ever need that information in your movement code for some reason.
- As I mentioned multiple times already, the iterative method for moving up close to a wall could be replaced by a slightly more sophisticated method that would directly find the farthest free position.


# Closing Words

I hope this tutorial made you more capable of solving these kinds of collision and movement challenges by yourself. You're free to use the code of this project for your own games.

## Further Reading
I recommend this 2012 article, [_The guide to implementing 2D platformers_](http://higherorderfun.com/blog/2012/05/20/the-guide-to-implementing-2d-platformers/) by Rodrigo Monteiro, which concisely explains how the simpler types of 2D platformers can be implemented, with an emphasis on tile-based platformers.

The [Sonic Physics Guide](https://info.sonicretro.org/Sonic_Physics_Guide) explains the mechanics of the original Sonic games in a lot of detail. It's a good resource, but I think the approach these games took is outdated. The amount of optimizations that today would offer insignificant gains has a huge cost on the flexibility of the engine.

--- 

Author: Pere Dolcet

Editor: Ted Bendixson

(02/2023)
