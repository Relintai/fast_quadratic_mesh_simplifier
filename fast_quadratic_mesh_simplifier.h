#ifndef FAST_QUADRATIC_MESH_SIMPLIFIER_H
#define FAST_QUADRATIC_MESH_SIMPLIFIER_H

/*

Copyright (c) 2020 Péter Magyar
Copyright(c) 2017-2020 Mattias Edlund

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "core/reference.h"

#include <limits>

#include "mesh_utils.h"

#include "core/array.h"
#include "core/pool_vector.h"
#include "core/resource.h"

class FastQuadraticMeshSimplifier : public Reference {
	GDCLASS(FastQuadraticMeshSimplifier, Reference);

public:
	void initialize(const Array &arrays);
	Array get_arrays();
	void simplify_mesh(float quality);
	void simplify_mesh_lossless();

	void update_mesh(int iteration);
	void refresh_vertices();
	void update_references();
	int remove_vertex_pass(int startTrisCount, int targetTrisCount, double threshold, PoolVector<bool> &deleted0, PoolVector<bool> &deleted1, int deletedTris);
	void compact_mesh();
	bool are_uvs_the_same(int channel, int indexA, int indexB);
	double vertex_error(SymmetricMatrix q, double x, double y, double z);
	double calculate_error(MUVertex vert0, MUVertex vert1, Vector3 *result);
	int update_triangles(int i0, int ia0, const MUVertex &v, PoolVector<bool> &deleted, int deletedTriangles);
	bool flipped(const Vector3 &p, int i0, int i1, const MUVertex &v0, PoolVector<bool> &deleted);
	static Vector3 calculate_barycentric_coords(Vector3 const &point, Vector3 const &a, Vector3 const &b, Vector3 const &c);
	void interpolate_vertex_attributes(int dst, int i0, int i1, int i2, Vector3 &barycentricCoord);

	static double min3(double val1, double val2, double val3) {
		return (val1 < val2 ? (val1 < val3 ? val1 : val3) : (val2 < val3 ? val2 : val3));
	}

	FastQuadraticMeshSimplifier();

protected:
	static void _bind_methods();

private:
	bool _has_normals;
	bool _has_colors;
	bool _has_uvs;
	bool _has_uv2s;
	PoolVector<int> _indices;

	PoolVector<MUTriangle> _mu_triangles;
	PoolVector<MUVertex> _mu_vertices;
	PoolVector<MURef> _mu_refs;

	double _vertex_link_distance_sqr = std::numeric_limits<double>::epsilon();
	int _max_iteration_count;
	double _agressiveness;
	bool _enable_smart_link;
	bool _preserve_border_dges;
	bool _preserve_uv_seam_edges;
	bool _preserve_uv_foldover_edges;
};

#endif