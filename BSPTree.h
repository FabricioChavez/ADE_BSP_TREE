#ifndef BSPTREE_H
#define BSPTREE_H

#include <memory>
#include <vector>
#include <functional>
#include "Plane.h"
#include "Ball.h"
#include <stack>
// Forward declarations
template <typename T>
class BSPNode;

template <typename T>
class BSPTree;

// BSPNode class template
template <typename T = NType>
class BSPNode {
private:
    Plane<T> partition_;
    std::unique_ptr<BSPNode<T>> front_;
    std::unique_ptr<BSPNode<T>> back_;
    std::vector<Polygon<T>> polygons_;
    static const T zero;

public:
    BSPNode() : partition_(), front_(nullptr), back_(nullptr) {}
    explicit BSPNode(const Polygon<T> & pol ): partition_(pol.getPlane()),front_(nullptr),back_(nullptr){
        polygons_.push_back(pol);
    }
    ~BSPNode() = default;

    BSPNode(const BSPNode&) = delete;
    BSPNode& operator=(const BSPNode&) = delete;

    // Getters
    const Plane<T>& getPartition() const { return partition_; }
    const std::vector<Polygon<T>>& getPolygons() const { return polygons_; }
    const BSPNode<T>* getFront() const { return front_.get(); }
    const BSPNode<T>*  getBack() const { return  back_.get(); }

    void insert(const Polygon<T>& polygon);

    // Método de consulta: recolecta en 'results' los polígonos que pueden colisionar con la Ball.
    void query(const Ball<T>& ball, const LineSegment<T>& movement, std::vector<Polygon<T>>& results) const;

    void candidate_tree(const Point3D<T>& pos, std::stack<const BSPNode<T>*>& path) const;

    void region_explorer(const BSPNode<T>* node, std::vector<Polygon<T>>& ans, const Ball<T>& ball, const LineSegment<T>& movement) const;

    bool is_there_interception(const Ball<T>& ball, const LineSegment<T>& movement, const Polygon<T>& poly) const;

    // Print
    void print(std::ostream& os, int indent = 0) const{
        std::string indentStr(indent * 4, ' ');

        os << indentStr << "BSPNode:\n";
        os << indentStr << "  Partition: " << partition_ << "\n";

        // Imprimir  polígonos
        os << indentStr << "  Polygons (" << polygons_.size() << "): ";
        if (polygons_.empty()) {
            os << "None\n";
        } else {
            os << "\n";
            for (size_t i = 0; i < polygons_.size(); ++i) {
                os << indentStr << "    [" << i << "]: " << polygons_[i] << "\n";
            }
        }

        // Front
        if (front_) {
            os << indentStr << "  Front:\n";
            front_->print(os, indent + 1);
        } else {
            os << indentStr << "  Front: NULL\n";
        }

        // Back
        if (back_) {
            os << indentStr << "  Back:\n";
            back_->print(os, indent + 1);
        } else {
            os << indentStr << "  Back: NULL\n";
        }
    }


    // Recorrido
    void collectNodes(std::vector<const BSPNode<T>*>& nodes) const{
        nodes.push_back(this);
        if (front_)
            front_->collectNodes(nodes);
        if (back_)
            back_->collectNodes(nodes);
    }

    void collectPolygons(std::vector<Polygon<T>>& polys) const {
        polys.insert(polys.end(), polygons_.begin(), polygons_.end());
        if (front_)
            front_->collectPolygons(polys);
        if (back_)
            back_->collectPolygons(polys);
    }

    void traverse(std::function<void(const BSPNode<T>&)> func) const {
        func(*this);
        if (front_)
            front_->traverse(func);
        if (back_)
            back_->traverse(func);
    }
};

template<typename T>
void BSPNode<T>::insert(const Polygon<T> &polygon) {

    RelationType relation = polygon.relationWithPlane(partition_);


    switch (relation) {

        case COINCIDENT :
            polygons_.push_back(polygon);
            break;

        case IN_FRONT:

            if(front_== nullptr){
                front_ = std::make_unique<BSPNode<T>>(polygon);
            }else front_->insert(polygon);

            break;
        case BEHIND:
            if(back_== nullptr){
                back_ = std::make_unique<BSPNode<T>>(polygon);
            }else back_->insert(polygon);

            break;
        case SPLIT:

            std::pair<Polygon<T>, Polygon<T>> ans = polygon.split(partition_);

            if(front_== nullptr){
                front_ = std::make_unique<BSPNode<T>>(ans.first);
            }else front_->insert(ans.first);

            if(back_== nullptr){
                back_ = std::make_unique<BSPNode<T>>(ans.second);
            }else back_->insert(ans.second);

            break;
    }
}

template <typename T>
const T BSPNode<T>::zero = static_cast<NType>(0.0);

template<typename T>
void BSPNode<T>::query(const Ball<T> &ball, const LineSegment<T> &movement, std::vector<Polygon<T>> &results) const {

    std::stack<const BSPNode<T>*> r0_path;
    std::stack<const BSPNode<T>*> rf_path;

    candidate_tree(movement.getP1()  , r0_path);
    candidate_tree(movement.getP2() , rf_path);


    while (r0_path.top()!=rf_path.top()){

        if(r0_path.size() == rf_path.size()){
            if(r0_path.top()!=rf_path.top()){
                r0_path.pop() ;
                rf_path.pop();
            }
        }else if(r0_path.size() > rf_path.size()) r0_path.pop();
        else if(r0_path.size() < rf_path.size()) rf_path.pop();
    }

    const BSPNode<T>* sub_tree_search_root = rf_path.top();




    region_explorer(sub_tree_search_root , results, ball , movement);



}


template<typename T>
void BSPNode<T>::candidate_tree(const Point3D<T>& pos, std::stack<const BSPNode<T>*>& path) const{
    if(back_ == nullptr and front_ == nullptr){
        path.push(this);
        return;
    }

    Plane<T> plane = getPartition();

    T loc = plane.distance(pos);

    if(loc== zero) {
        path.push(this);
        return;
    }

    path.push(this);

    if(loc > zero and front_!= nullptr) {
        front_->candidate_tree(pos , path);
    }else if(loc < zero and back_ != nullptr){
        back_->candidate_tree(pos , path);
    }
}

template<typename T>
void BSPNode<T>::region_explorer(const BSPNode<T>* node, std::vector<Polygon<T>>& ans, const Ball<T>& ball, const LineSegment<T>& movement) const {

    for(auto p : node->polygons_){
        if(is_there_interception(ball, movement , p))
            ans.push_back(p);
    }

    if(node->getBack()== nullptr and node->getFront() == nullptr) return;

    if(node->front_ != nullptr) region_explorer(node->getFront() , ans , ball , movement);
    if(node->back_ != nullptr) region_explorer(node->getBack() , ans , ball , movement);

}

template<typename T>
bool BSPNode<T>::is_there_interception(const Ball<T> &ball, const LineSegment<T> &movement, const Polygon<T> &poly) const {
    Plane<T> plane = poly.getPlane();
    T r = ball.getRadius();
    T dStart = plane.distance(movement.getP1());
    T dEnd   = plane.distance(movement.getP2());

    // Claro que no hay interseccion
    if ((dStart > r && dEnd > r) || (dStart < -r && dEnd < -r))
        return false;

    // Calcular t (si existe)
    T denom = dStart - dEnd;
    if (denom == 0) return false; // Evitar división por cero.
    T t = dStart / denom;
    if (t < T(0) || t > T(1))
        return false;

    // Punto de interseccion
    Point3D<T> intersection = movement.getP1() + (movement.getP2() - movement.getP1()) * t;
    return poly.contains(intersection);
}

// BSPTree class template
template <typename T = NType>
class BSPTree {
private:
    std::unique_ptr<BSPNode<T>> root_;

public:
    BSPTree() : root_(nullptr) {}
    ~BSPTree() = default;

    void insert(const Polygon<T>& polygon);


    // Devuelve los polígonos candidatos a colisión con la Ball.
    std::vector<Polygon<T>> query(const Ball<T>& ball, const LineSegment<T>& movement) const;

    // Print
    void print(std::ostream& os) const{
        if (root_) {
            os << "BSPTree:\n";
            root_->print(os, 1);
        } else {
            os << "BSPTree is empty.\n";
        }
    }


    // Funciones de recorrido
    std::vector<const BSPNode<T>*> getAllNodes() const{
        std::vector<const BSPNode<T>*> nodes;
        if (root_)
            root_->collectNodes(nodes);
        return nodes;
    }
    std::vector<Polygon<T>> getAllPolygons() const{
        std::vector<Polygon<T>> polys;
        if (root_)
            root_->collectPolygons(polys);
        return polys;
    }
    void traverse(std::function<void(const BSPNode<T>&)> func) const{
        if (root_)
            root_->traverse(func);
    }
};

template<typename T>
void BSPTree<T>::insert(const Polygon<T> &polygon) {
    if(root_ == nullptr){
        root_ = std::make_unique<BSPNode<T>>(polygon);
    }else {
        root_->insert(polygon);
    }
}

template<typename T>
std::vector<Polygon<T>> BSPTree<T>::query(const Ball<T> &ball, const LineSegment<T> &movement) const {
    std::vector<Polygon<T>> ans;

    root_->query(ball , movement , ans);
    return ans;
}


#endif // BSPTREE_H
