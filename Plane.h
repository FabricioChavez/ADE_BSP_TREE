#ifndef PLANE_H
#define PLANE_H

#include "DataType.h"
#include "Point.h"
#include "Line.h"
#include <vector>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <cmath>
using namespace std;

// Forward declarations
template <typename T>
class Plane;

template <typename T>
class Polygon;

// Plane class template
template <typename T = NType>
class Plane {
private:
    Point3D <T>  point_;    // A point on the plane
    Vector3D<T> normal_;    // A normal vector to the plane
    T const zero_scalar ;
    Vector3D<T> const zero_vector;
public:
    // Constructors
    Plane() : point_(),
              normal_(Vector3D<T>(static_cast<T>(0), static_cast<T>(0), static_cast<T>(1))),
              zero_scalar(static_cast<NType>(0.0)),
              zero_vector(static_cast<NType>(0.0), static_cast<NType>(0.0), static_cast<NType>(0.0)){

    }
    Plane(const Point3D<T>& point, const Vector3D<T>& normal) : point_(point), normal_(normal.normalized()),
                                                                zero_scalar(static_cast<T>(0.0)),
                                                                zero_vector(static_cast<T>(0.0), static_cast<T>(0.0), static_cast<T>(0.0)){

    }

    // Distance from a point to the plane
    T distance(const Point3D<T>& p) const ;

    // Intersection with a line
    Point3D<T> intersect(const Line<T>& l) const;

    // Containment checks
    bool contains(const Point3D<T>& p) const;
    bool contains(const Line<T>& l) const;

    // Getters
    Point3D<T> getPoint() const { return point_; }
    Vector3D<T> getNormal() const { return normal_; }

    // Setters
    void setPoint(const Point3D<T>& point) { point_ = point; }
    void setNormal(const Vector3D<T>& normal) { normal_ = normal.normalized(); }

    // Operators
    bool operator==(const Plane& other) const;
    bool operator!=(const Plane& other) const { return !(*this == other); }

    // Output operator
    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Plane<U>& plane);
};

// Polygon class template
template <typename T = NType>
class Polygon {
private:
    std::vector<Point3D<T>> vertices_;
    static const T  zero_scalar ;
    static const Vector3D<T>  zero_vector;

public:
    // Constructors
    Polygon() : vertices_() {}
    Polygon(const std::vector<Point3D<T>>& vertices) : vertices_(vertices) {}

    // Getters
    std::vector<Point3D<T>> getVertices() const { return vertices_; }
    const Point3D<T>& getVertex(size_t index) const { return vertices_.at(index); }
    Plane<T>    getPlane   () const;    // Get the plane    of the polygon
    Vector3D<T> getNormal  () const;    // Get the normal   of the polygon
    Point3D<T>  getCentroid() const;    // Get the centroid of the polygon
    Point3D<T>  getRandomPoint() ; // Get a random point contained in the polygon plane

    // Setters
    void setVertices(const std::vector<Point3D<T>>& vertices) { vertices_ = vertices; }

    // Check if a point is inside the polygon
    bool contains(const Point3D<T>& p) const;

    // Get the relation of the polygon with a plane
    RelationType relationWithPlane(const Plane<T>& plane) const;

    // Split the polygon by a plane
    std::pair<Polygon<T>, Polygon<T>> split(const Plane<T>& plane) const;

    // Compute the area of the polygon
    T area() const;

    // Operators
    bool operator==(const Polygon& other) const;
    bool operator!=(const Polygon& other) const { return !(*this == other); }



    // Output operator
    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Polygon<U>& polygon);
};

//Constants
template <typename T>
const T Polygon<T>::zero_scalar = static_cast<NType>(0.0);

template <typename T>
const Vector3D<T> Polygon<T>::zero_vector = Vector3D<T>(static_cast<NType>(0.0), static_cast<NType>(0.0), static_cast<NType>(0.0));

// Equality operators
template <typename T>
bool Plane<T>::operator==(const Plane<T>& other) const {
    bool normalsEqual = (normal_ == other.normal_) || (normal_ == -other.normal_);
    return normalsEqual && contains(other.point_);
}

// Output operator for Plane
template <typename T>
std::ostream& operator<<(std::ostream& os, const Plane<T>& plane) {
    os << "Point: " << plane.point_ << ", Normal: " << plane.normal_;
    return os;
}

template<typename T>
T Plane<T>::distance(const Point3D<T> &p) const { // OK DONE
    Vector3D <T> plane_to_point(p - point_);
    normal_.normalized();
    return plane_to_point.dotProduct(normal_);
}

template<typename T>
Point3D<T> Plane<T>::intersect(const Line<T> &l) const { // ok DONE

    // if(contains(l)) throw std::runtime_error("the line is on the plane");

    //  if(l.isOrthogonal(getNormal()) and distance(l.getPoint()) > zero_scalar) throw std::runtime_error("the line does not have any interception");

    Point3D<T> plane_to_point = point_ - l.getPoint();


    T intersect_parameter  =  normal_.dotProduct(plane_to_point)/normal_.dotProduct(l.getDirection());

    Vector3D<T> m = l.getDirection();

    //  std::cout<<"INTERSECT PARAMETER "<<intersect_parameter<<std::endl;

    Point3D<T> movement(m.getX() , m.getY() , m.getZ());
    movement*=intersect_parameter;

    //  std::cout<<" MOVEMENT "<<movement<<std::endl;


    return l.getPoint() + movement ; // line vectorial equation
}



template<typename T>
bool Plane<T>::contains(const Point3D<T> &p) const { // ok done
    return distance(p)==zero_scalar;
}

template<typename T>
bool Plane<T>::contains(const Line<T> &l) const { // OK DONE



    T dist = distance(l.getPoint());
    std::cout<<"DISTANCE IS "<< dist<<std::endl;
    std::cout<<"IS PARALLEL TO PLANE ? " << std::boolalpha<< l.isOrthogonal(getNormal())<< std::endl;

    return l.isOrthogonal(getNormal()) and  dist == zero_scalar  ;


}

// Equality operators
template <typename T>
bool Polygon<T>::operator==(const Polygon<T>& other) const {
    return vertices_ == other.vertices_;
}

// Output operator for Polygon
template <typename T>
std::ostream& operator<<(std::ostream& os, const Polygon<T>& polygon) {
    os << "Vertices: ";
    for (const auto& vertex : polygon.vertices_) {
        os << vertex << " ";
    }
    return os;
}



template<typename T>
Plane<T> Polygon<T>::getPlane() const { // ok done

    Point3D<T> origin = vertices_[0];
    Point3D<T> pa = vertices_[1];
    Point3D<T> pb = vertices_[2];

    Vector3D<T> v_a(pa - origin);
    Vector3D<T> v_b(pb - origin);
    Vector3D<T> normal = v_a.crossProduct(v_b);

    return  Plane<T> (origin , normal);


}

template<typename T>
Vector3D<T> Polygon<T>::getNormal() const {
    Point3D<T> origin = vertices_[0];
    Point3D<T> pa = vertices_[1];
    Point3D<T> pb = vertices_[2];

    Vector3D<T> v_a(pa - origin);
    Vector3D<T> v_b(pb - origin);
    Vector3D<T> normal = v_a.crossProduct(v_b);
    normal.normalize();

    return  normal;
}

template<typename T>
Point3D<T> Polygon<T>::getCentroid() const { //ok done
    Point3D<T> zero_point(0.0,0.0,0.0);

    for(auto p : vertices_){
        zero_point+=p;
    }
    zero_point/=(vertices_.size());
    return  zero_point;
}

template<typename T>
bool Polygon<T>::contains(const Point3D<T> &p) const { // ok done

    std::vector<Point3D<T>> vertices = vertices_;

    vertices.push_back(vertices_[0]);
    Vector3D<T> v1(vertices[0]-p);
    Vector3D<T> v2(vertices[1]-p);

    Vector3D<T> cruz = v2.crossProduct(v1);
    T d = cruz.dotProduct(getNormal());
    //   cout<<"NORMAL "<<getNormal()<<endl;
    bool sign = d > zero_scalar;
    //  cout<<"SIGN FIRST "<<sign<<endl;
    int n = int(vertices.size());
    for (int i = 2; i < n; ++i) {
        v1= Vector3D<T>(vertices[i-1]-p);
        v2 = Vector3D<T>(vertices[i]-p);
        //cout<<vertices[0]<< "PREVIUS"<<endl;
        //cout<<v2<< "CURRENT"<<endl;
        cruz = v2.crossProduct(v1);
        // cout<<"NEW NORMAL "<<cruz<<endl;
        d = cruz.dotProduct(getNormal());
        // cout<<d<<endl;
        bool new_sign = d > zero_scalar;
        //cout<<"SIGN --->"<<new_sign<<endl;
        if(new_sign != sign) return false;

    }
    return true;

}

template<typename T>
RelationType Polygon<T>::relationWithPlane(const Plane<T> &plane) const {

    int n = getVertices().size();

    int pos_vertice = 0;
    int neg_vertice = 0;

    for (int i = 0; i < n; ++i) {
        Vector3D<T> vertice_to_plane(getVertices()[i]- plane.getPoint());

        T location = vertice_to_plane.dotProduct(plane.getNormal());

        if(location > zero_scalar) pos_vertice++;
        else if(location < zero_scalar) neg_vertice++;
    }

    if(pos_vertice == 0 and neg_vertice == 0) return RelationType::COINCIDENT;


    if(pos_vertice > 0 and neg_vertice > 0) return RelationType::SPLIT;

    if(pos_vertice> 0 and neg_vertice == 0) return  RelationType::IN_FRONT;

    if(pos_vertice == 0  and neg_vertice > 0) return RelationType::BEHIND;

   // std::cout<<"It should not reach this place , relationship not defined "<<std::endl;
    return RelationType::COINCIDENT;

}
template<typename T>
std::pair<Polygon<T>, Polygon<T>> Polygon<T>::split(const Plane<T> &plane) const { // ok done

    std::vector<Point3D<T>> vertices = vertices_;;//push_back to cover all edges in the polygon
    vertices.push_back(vertices_[0]);
    std::vector<Point3D<T>> pos_p , neg_p;

    int n = vertices.size();
    for (int i = 1; i < n; ++i) {
        Point3D<T> vertice_to_plane_c(vertices[i]); //current point location
        Point3D<T> vertice_to_plane_b(vertices[i-1]); //before point location

        //  cout<<"CURRENT "<< vertices[i]<<endl;
        //  cout<<"BEFORE  "<<vertices[i-1]<<endl;

        T location_c = plane.distance(vertice_to_plane_c);
        T location_b = plane.distance(vertice_to_plane_b);

        T dist= abs(plane.distance(vertices[i]));
        T dist_b= abs(plane.distance(vertices[i-1]));
        if((dist==zero_scalar))
        {
            pos_p.push_back(vertices[i]);
            neg_p.push_back(vertices[i]);
            continue;
            // cout<<"CURRENT IS ON THE PLANE -_-"<<boolalpha<<(dist==zero_scalar)<<endl;
        }




        // front 1  back 0
        bool current_pos_c = (location_c > zero_scalar);
        bool current_pos_b = (location_b > zero_scalar);





        if(current_pos_c ^ current_pos_b and dist_b > zero_scalar){ // if both are on different sides of the plane there is an intersection in this segment

            Line<T> inter_seg(vertices[i-1],vertices[i]);
            Point3D<T> intersect_point = plane.intersect(inter_seg);
            //   cout<<"INSERSECT"<<intersect_point<<endl;

            neg_p.push_back(intersect_point);
            pos_p.push_back(intersect_point);
            if(current_pos_c ){
                pos_p.push_back(vertices[i]);
            }else {
                neg_p.push_back(vertices[i]);

            }


            continue;

        }


        if(current_pos_c) pos_p.push_back(vertices[i]);
        else neg_p.push_back(vertices[i]);




    }

    Polygon<T> pos_poly(pos_p);
    Polygon<T> neg_poly(neg_p);



    return std::make_pair(pos_poly , neg_poly);
}
template<typename T>
T Polygon<T>::area() const { // ok done

    int n = vertices_.size();
    Point3D<T> origin = vertices_[0];

    T area = static_cast<T> (0.0);

    Point3D<T> init(0,0,0);

    for (int i = 2; i < n; ++i) {
        Vector3D<T> va(vertices_[i]-origin);
        Vector3D<T> vb(vertices_[i-1]-origin);

        Vector3D<T> cross = va.crossProduct(vb);

        Point3D<T> add(cross.getX() , cross.getY() , cross.getZ());
        init+=add;


    }

    area = init.magnitude();


    return area*static_cast<T>(0.5);
}

template<typename T>
Point3D<T> Polygon<T>::getRandomPoint() {

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<float> centerDist(-100.0f, 100.0f);
    T alfa = static_cast<T>(centerDist(gen));
    T beta = static_cast<T>(centerDist(gen));

    Vector3D<T> u(vertices_[0], vertices_[1]);
    Vector3D<T> v(vertices_[0], vertices_[2]);

    Point3D<T> u_(u.getX(), u.getY() , u.getZ());
    Point3D<T> v_(v.getX(), v.getY() , v.getZ());

    return  vertices_[0]+alfa*u_+beta*v_;

}

#endif // PLANE_H