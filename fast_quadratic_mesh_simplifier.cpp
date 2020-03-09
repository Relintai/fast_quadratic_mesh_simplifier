#include "fast_quadratic_mesh_simplifier.h"

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

#include "scene/resources/mesh.h"
#include "servers/visual_server.h"

int FastQuadraticMeshSimplifier::get_max_iteration_count() const {
	return _max_iteration_count;
}
void FastQuadraticMeshSimplifier::set_max_iteration_count(const int value) {
	_max_iteration_count = value;
}

double FastQuadraticMeshSimplifier::get_agressiveness() const {
	return _agressiveness;
}
void FastQuadraticMeshSimplifier::set_agressiveness(const double value) {
	_agressiveness = value;
}

bool FastQuadraticMeshSimplifier::get_enable_smart_link() const {
	return _enable_smart_link;
}
void FastQuadraticMeshSimplifier::set_enable_smart_link(const bool value) {
	_enable_smart_link = value;
}

bool FastQuadraticMeshSimplifier::get_preserve_border_dges() const {
	return _preserve_border_dges;
}
void FastQuadraticMeshSimplifier::set_preserve_border_dges(const bool value) {
	_preserve_border_dges = value;
}

bool FastQuadraticMeshSimplifier::get_preserve_uv_seam_edges() const {
	return _preserve_uv_seam_edges;
}
void FastQuadraticMeshSimplifier::set_preserve_uv_seam_edges(const bool value) {
	_preserve_uv_seam_edges = value;
}

bool FastQuadraticMeshSimplifier::get_preserve_uv_foldover_edges() const {
	return _preserve_uv_foldover_edges;
}
void FastQuadraticMeshSimplifier::set_preserve_uv_foldover_edges(const bool value) {
	_preserve_uv_foldover_edges = value;
}

void FastQuadraticMeshSimplifier::initialize(const Array &arrays) {
	ERR_FAIL_COND(arrays.size() != ArrayMesh::ARRAY_MAX);

	PoolVector<Vector3> vertices = arrays.get(ArrayMesh::ARRAY_VERTEX);
	PoolVector<Vector3> normals = arrays.get(ArrayMesh::ARRAY_NORMAL);
	PoolVector<Color> colors = arrays.get(ArrayMesh::ARRAY_COLOR);
	PoolVector<Vector2> uvs = arrays.get(ArrayMesh::ARRAY_TEX_UV);
	PoolVector<Vector2> uv2s = arrays.get(ArrayMesh::ARRAY_TEX_UV2);

	_format = 0;

	if (normals.size() > 0)
		_format |= VisualServer::ARRAY_FORMAT_NORMAL;

	if (colors.size() > 0)
		_format |= VisualServer::ARRAY_FORMAT_COLOR;

	if (uvs.size() > 0)
		_format |= VisualServer::ARRAY_FORMAT_TEX_UV;

	if (uv2s.size() > 0)
		_format |= VisualServer::ARRAY_FORMAT_TEX_UV2;

	_vertices.resize(vertices.size());
	for (int i = 0; i < vertices.size(); ++i) {

		Vertex vert;
		vert.vertex = vertices[i];

		if (normals.size() > i)
			vert.normal = normals[i];

		if (colors.size() > i)
			vert.color = colors[i];

		if (uvs.size() > i)
			vert.uv = uvs[i];

		if (uv2s.size() > i)
			vert.uv2 = uv2s[i];

		_vertices.set(i, vert);
	}

	PoolVector<int> indices = arrays.get(ArrayMesh::ARRAY_INDEX);
	_indices.resize(indices.size());
	for (int i = 0; i < indices.size(); ++i) {
		_indices.set(i, indices[i]);
	}

	if ((_indices.size() % 3) != 0)
		ERR_FAIL_MSG("The index array length must be a multiple of 3 in order to represent triangles.");

	int triangle_count = _indices.size() / 3;
	_mu_triangles.resize(triangle_count);

	for (int i = 0; i < triangle_count; ++i) {
		int offset = i * 3;
		int v0 = _indices[offset];
		int v1 = _indices[offset + 1];
		int v2 = _indices[offset + 2];
		_mu_triangles.set(i, MUTriangle(v0, v1, v2, 0));
	}

	_mu_vertices.resize(_vertices.size());
	for (int i = 0; i < _vertices.size(); ++i) {
		_mu_vertices.set(i, MUVertex(_vertices[i]));
	}
}

Array FastQuadraticMeshSimplifier::get_arrays() {
	Array arr;

	arr.resize(ArrayMesh::ARRAY_MAX);

	PoolVector<Vector3> vertices;
	PoolVector<Vector3> normals;
	PoolVector<Color> colors;
	PoolVector<Vector2> uvs;
	PoolVector<Vector2> uv2s;
	PoolVector<int> indices;

	vertices.resize(_mu_vertices.size());

	if ((_format & VisualServer::ARRAY_FORMAT_NORMAL) != 0)
		normals.resize(_mu_vertices.size());

	if ((_format & VisualServer::ARRAY_FORMAT_COLOR) != 0)
		colors.resize(_mu_vertices.size());

	if ((_format & VisualServer::ARRAY_FORMAT_TEX_UV) != 0)
		uvs.resize(_mu_vertices.size());

	if ((_format & VisualServer::ARRAY_FORMAT_TEX_UV2) != 0)
		uv2s.resize(_mu_vertices.size());

	indices.resize(_indices.size());

	for (int i = 0; i < vertices.size(); ++i) {
		vertices.set(i, _mu_vertices[i].vertex.vertex);
	}

	for (int i = 0; i < normals.size(); ++i) {
		normals.set(i, _mu_vertices[i].vertex.normal);
	}

	for (int i = 0; i < colors.size(); ++i) {
		colors.set(i, _mu_vertices[i].vertex.color);
	}

	for (int i = 0; i < uvs.size(); ++i) {
		uvs.set(i, _mu_vertices[i].vertex.uv);
	}

	for (int i = 0; i < uv2s.size(); ++i) {
		uv2s.set(i, _mu_vertices[i].vertex.uv2);
	}

	for (int i = 0; i < indices.size(); ++i) {
		indices.set(i, _indices[i]);
	}

	arr.set(ArrayMesh::ARRAY_VERTEX, vertices);
	arr.set(ArrayMesh::ARRAY_NORMAL, normals);
	arr.set(ArrayMesh::ARRAY_COLOR, colors);
	arr.set(ArrayMesh::ARRAY_TEX_UV, uvs);
	arr.set(ArrayMesh::ARRAY_TEX_UV2, uv2s);
	arr.set(ArrayMesh::ARRAY_INDEX, indices);

	return arr;
}

void FastQuadraticMeshSimplifier::simplify_mesh(float quality) {
	quality = CLAMP(quality, 0, 1);

	int deletedTris = 0;
	PoolVector<bool> deleted0;
	deleted0.resize(20);
	PoolVector<bool> deleted1;
	deleted1.resize(20);

	int startTrisCount = _mu_triangles.size();
	int targetTrisCount = static_cast<int>(_mu_triangles.size() * quality + 0.5);

	for (int iteration = 0; iteration < _max_iteration_count; iteration++) {
		if ((startTrisCount - deletedTris) <= targetTrisCount)
			break;

		// Update mesh once in a while
		if ((iteration % 5) == 0) {
			update_mesh(iteration);
		}

		// Clear dirty flag
		for (int i = 0; i < _mu_triangles.size(); ++i) {
			MUTriangle t = _mu_triangles[i];

			t.dirty = false;

			_mu_triangles.set(i, t);
		}

		// All triangles with edges below the threshold will be removed
		//
		// The following numbers works well for most models.
		// If it does not, try to adjust the 3 parameters
		double threshold = 0.000000001 * Math::pow(iteration + 3, _agressiveness);

		print_error("iteration " + String::num(iteration) + " - triangles " + String::num((startTrisCount - deletedTris)) + " threshold " + String::num(threshold));

		// Remove vertices & mark deleted triangles
		remove_vertex_pass(startTrisCount, targetTrisCount, threshold, &deleted0, &deleted1, &deletedTris);
	}

	compact_mesh();

	print_error("Finished simplification with triangle count " + String::num(_mu_triangles.size()));
}

void FastQuadraticMeshSimplifier::simplify_mesh_lossless() {
	int deletedTris = 0;
	PoolVector<bool> deleted0;
	PoolVector<bool> deleted1;
	int startTrisCount = _mu_triangles.size();

	for (int iteration = 0; iteration < 9999; iteration++) {
		// Update mesh constantly
		update_mesh(iteration);

		// Clear dirty flag
		for (int i = 0; i < _mu_triangles.size(); ++i) {
			MUTriangle t = _mu_triangles[i];

			t.dirty = false;

			_mu_triangles.set(i, t);
		}

		// All triangles with edges below the threshold will be removed
		//
		// The following numbers works well for most models.
		// If it does not, try to adjust the 3 parameters
		double threshold = 1.0E-3;

		//Debug.LogFormat("Lossless iteration {0} - triangles {1}", iteration, triangleCount);

		// Remove vertices & mark deleted triangles
		remove_vertex_pass(startTrisCount, 0, threshold, &deleted0, &deleted1, &deletedTris);

		if (deletedTris <= 0)
			break;

		deletedTris = 0;
	}

	compact_mesh();

	//Debug.LogFormat("Finished simplification with triangle count {0}", this.triangles.Length);
}

void FastQuadraticMeshSimplifier::update_mesh(int iteration) {
	if (iteration > 0) // compact triangles
	{
		int dst = 0;
		for (int i = 0; i < _mu_triangles.size(); ++i) {
			MUTriangle t = _mu_triangles[i];

			if (!t.deleted) {
				if (dst != i) {
					_mu_triangles.set(dst, t);
				}

				++dst;
			}
		}

		_mu_triangles.resize(dst);
	}

	update_references();

	// Identify boundary : vertices[].border=0,1
	if (iteration == 0) {
		PoolVector<int> vcount;
		PoolVector<int> vids;

		int vsize = 0;
		for (int i = 0; i < _mu_vertices.size(); ++i) {
			MUVertex v = _mu_vertices[i];

			v.border_edge = false;
			v.uv_seam_edge = false;
			v.uv_foldover_edge = false;

			_mu_vertices.set(i, v);
		}

		int ofs = 0;
		int id = 0;
		int borderVertexCount = 0;
		double borderMinX = std::numeric_limits<double>::max();
		double borderMaxX = std::numeric_limits<double>::min();

		for (int i = 0; i < _mu_vertices.size(); i++) {
			MUVertex v = _mu_vertices[i];

			int tstart = v.tstart;
			int tcount = v.tcount;
			vcount.resize(0);
			vids.resize(0);
			vsize = 0;

			for (int j = 0; j < tcount; ++j) {
				int tid = _mu_refs[tstart + j].tid;
				MUTriangle t = _mu_triangles[tid];

				for (int k = 0; k < 3; ++k) {
					ofs = 0;
					id = t.get(k);

					while (ofs < vsize) {
						if (vids[ofs] == id)
							break;

						++ofs;
					}

					if (ofs == vsize) {
						vcount.push_back(1);
						vids.push_back(id);
						++vsize;
					} else {
						vcount.set(ofs, vcount[ofs] + 1);
					}
				}
			}

			for (int j = 0; j < vsize; ++j) {
				if (vcount[j] == 1) {
					id = vids[j];

					MUVertex vv = _mu_vertices[id];
					vv.border_edge = true;
					_mu_vertices.set(id, vv);

					++borderVertexCount;

					if (_enable_smart_link) {
						if (v.vertex.vertex.x < borderMinX) {
							borderMinX = v.vertex.vertex.x;
						}

						if (v.vertex.vertex.x > borderMaxX) {
							borderMaxX = v.vertex.vertex.x;
						}
					}
				}
			}
		}

		if (_enable_smart_link) {
			// First find all border vertices
			Vector<BorderVertex> borderVertices;
			borderVertices.resize(borderVertexCount);

			int borderIndexCount = 0;
			double borderAreaWidth = borderMaxX - borderMinX;

			for (int i = 0; i < _mu_vertices.size(); ++i) {
				MUVertex mvi = _mu_vertices[i];

				if (mvi.border_edge) {
					int vertexHash = (int)(((((mvi.vertex.vertex.x - borderMinX) / borderAreaWidth) * 2.0) - 1.0) * std::numeric_limits<int>::max());
					borderVertices.set(borderIndexCount, BorderVertex(i, vertexHash));
					++borderIndexCount;
				}
			}

			// Sort the border vertices by hash
			borderVertices.sort_custom<BorderVertexComparer>();

			// Calculate the maximum hash distance based on the maximum vertex link distance
			double vertexLinkDistance = Math::sqrt(_vertex_link_distance_sqr);
			int hashMaxDistance = MAX((int)((vertexLinkDistance / borderAreaWidth) * std::numeric_limits<int>::max()), 1);

			// Then find identical border vertices and bind them together as one
			for (int i = 0; i < borderIndexCount; ++i) {
				int myIndex = borderVertices[i].index;

				if (myIndex == -1)
					continue;

				Vector3 myPoint = _mu_vertices[myIndex].vertex.vertex;

				for (int j = i + 1; j < borderIndexCount; j++) {

					int otherIndex = borderVertices[j].index;
					if (otherIndex == -1)
						continue;
					else if ((borderVertices[j].hash - borderVertices[i].hash) > hashMaxDistance) // There is no point to continue beyond this point
						break;

					Vector3 otherPoint = _mu_vertices[otherIndex].vertex.vertex;
					double sqrX = ((myPoint.x - otherPoint.x) * (myPoint.x - otherPoint.x));
					double sqrY = ((myPoint.y - otherPoint.y) * (myPoint.y - otherPoint.y));
					double sqrZ = ((myPoint.z - otherPoint.z) * (myPoint.z - otherPoint.z));
					double sqrMagnitude = sqrX + sqrY + sqrZ;

					if (sqrMagnitude <= _vertex_link_distance_sqr) {
						BorderVertex b = borderVertices[j]; // NOTE: This makes sure that the "other" vertex is not processed again
						b.index = -1;
						borderVertices.set(j, b);

						MUVertex miv = _mu_vertices[myIndex];
						MUVertex oiv = _mu_vertices[otherIndex];

						miv.border_edge = false;
						oiv.border_edge = false;

						if (are_uvs_the_same(0, myIndex, otherIndex)) {
							miv.uv_foldover_edge = true;
							oiv.uv_foldover_edge = true;
						} else {
							miv.uv_seam_edge = true;
							oiv.uv_seam_edge = true;
						}

						_mu_vertices.set(myIndex, miv);
						_mu_vertices.set(otherIndex, oiv);

						int otherTriangleCount = oiv.tcount;
						int otherTriangleStart = oiv.tstart;

						for (int k = 0; k < otherTriangleCount; k++) {
							MURef r = _mu_refs[otherTriangleStart + k];
							MUTriangle t = _mu_triangles[r.tid];
							t.set(r.tvertex, myIndex);

							_mu_triangles.set(r.tid, t);
						}
					}
				}
			}

			// Update the references again
			update_references();
		}

		// Init Quadrics by Plane & Edge Errors
		//
		// required at the beginning ( iteration == 0 )
		// recomputing during the simplification is not required,
		// but mostly improves the result for closed meshes
		for (int i = 0; i < _mu_vertices.size(); ++i) {
			MUVertex v = _mu_vertices[i];

			v.q.reset();

			_mu_vertices.set(i, v);
		}

		int v0, v1, v2;
		Vector3 n, p0, p1, p2, p10, p20, dummy;
		SymmetricMatrix sm;
		for (int i = 0; i < _mu_triangles.size(); ++i) {
			MUTriangle t = _mu_triangles[i];

			v0 = t.v0;
			v1 = t.v1;
			v2 = t.v2;

			MUVertex vt0 = _mu_vertices[v0];
			MUVertex vt1 = _mu_vertices[v1];
			MUVertex vt2 = _mu_vertices[v2];

			p0 = vt0.vertex.vertex;
			p1 = vt1.vertex.vertex;
			p2 = vt2.vertex.vertex;
			p10 = p1 - p0;
			p20 = p2 - p0;

			n = p10.cross(p20);

			n.normalize();

			t.n = n;

			_mu_triangles.set(i, t);

			sm.from_plane(n.x, n.y, n.z, -n.dot(p0));

			vt0.q += sm;
			vt1.q += sm;
			vt2.q += sm;

			_mu_vertices.set(v0, vt0);
			_mu_vertices.set(v1, vt1);
			_mu_vertices.set(v2, vt2);
		}

		for (int i = 0; i < _mu_triangles.size(); ++i) {
			// Calc Edge Error
			MUTriangle triangle = _mu_triangles[i];
			triangle.err0 = calculate_error(_mu_vertices[triangle.v0], _mu_vertices[triangle.v1], &dummy);
			triangle.err1 = calculate_error(_mu_vertices[triangle.v1], _mu_vertices[triangle.v2], &dummy);
			triangle.err2 = calculate_error(_mu_vertices[triangle.v2], _mu_vertices[triangle.v0], &dummy);
			triangle.err3 = FastQuadraticMeshSimplifier::min3(triangle.err0, triangle.err1, triangle.err2);
			_mu_triangles.set(i, triangle);
		}
	}
}

void FastQuadraticMeshSimplifier::update_references() {
	// Init Reference ID list
	for (int i = 0; i < _mu_vertices.size(); i++) {
		MUVertex v = _mu_vertices[i];

		v.tstart = 0;
		v.tcount = 0;

		_mu_vertices.set(i, v);
	}

	for (int i = 0; i < _mu_triangles.size(); i++) {

		MUTriangle t = _mu_triangles[i];

		MUVertex v = _mu_vertices[t.v0];
		++v.tcount;
		_mu_vertices.set(t.v0, v);

		v = _mu_vertices[t.v1];
		++v.tcount;
		_mu_vertices.set(t.v1, v);

		v = _mu_vertices[t.v2];
		++v.tcount;
		_mu_vertices.set(t.v2, v);
	}

	int tstart = 0;
	for (int i = 0; i < _mu_vertices.size(); ++i) {
		MUVertex v = _mu_vertices[i];

		v.tstart = tstart;
		tstart += v.tcount;
		v.tcount = 0;
		_mu_vertices.set(i, v);
	}

	// Write References
	_mu_refs.resize(tstart);
	for (int i = 0; i < _mu_triangles.size(); ++i) {
		MUTriangle t = _mu_triangles[i];

		int start0 = _mu_vertices[t.v0].tstart;
		int count0 = _mu_vertices[t.v0].tcount;
		int start1 = _mu_vertices[t.v1].tstart;
		int count1 = _mu_vertices[t.v1].tcount;
		int start2 = _mu_vertices[t.v2].tstart;
		int count2 = _mu_vertices[t.v2].tcount;

		MURef ref = _mu_refs[start0 + count0];
		ref.Set(i, 0);
		_mu_refs.set(start0 + count0, ref);

		ref = _mu_refs[start1 + count1];
		ref.Set(i, 1);
		_mu_refs.set(start1 + count1, ref);

		ref = _mu_refs[start2 + count2];
		ref.Set(i, 2);
		_mu_refs.set(start2 + count2, ref);

		MUVertex v = _mu_vertices[t.v0];
		++v.tcount;
		_mu_vertices.set(t.v0, v);

		v = _mu_vertices[t.v1];
		++v.tcount;
		_mu_vertices.set(t.v1, v);

		v = _mu_vertices[t.v2];
		++v.tcount;
		_mu_vertices.set(t.v2, v);
	}
}

/// <summary>
/// Finally compact mesh before exiting.
/// </summary>
void FastQuadraticMeshSimplifier::compact_mesh() {
	int dst = 0;

	for (int i = 0; i < _mu_vertices.size(); ++i) {
		MUVertex v = _mu_vertices.get(i);
		v.tcount = 0;
		_mu_vertices.set(i, v);
	}

	for (int i = 0; i < _mu_triangles.size(); i++) {
		MUTriangle triangle = _mu_triangles[i];

		if (!triangle.deleted) {
			if (triangle.va0 != triangle.v0) {
				int iDest = triangle.va0;
				int iSrc = triangle.v0;

				MUVertex d = _mu_vertices[iDest];
				MUVertex s = _mu_vertices[iSrc];
				d.vertex = s.vertex;
				_mu_vertices.set(iDest, d);

				triangle.v0 = triangle.va0;
			}

			if (triangle.va1 != triangle.v1) {
				int iDest = triangle.va1;
				int iSrc = triangle.v1;

				MUVertex d = _mu_vertices[iDest];
				MUVertex s = _mu_vertices[iSrc];
				d.vertex = s.vertex;
				_mu_vertices.set(iDest, d);

				triangle.v1 = triangle.va1;
			}

			if (triangle.va2 != triangle.v2) {
				int iDest = triangle.va2;
				int iSrc = triangle.v2;

				MUVertex d = _mu_vertices[iDest];
				MUVertex s = _mu_vertices[iSrc];
				d.vertex = s.vertex;
				_mu_vertices.set(iDest, d);

				triangle.v2 = triangle.va2;
			}

			int newTriangleIndex = ++dst;
			_mu_triangles.set(newTriangleIndex, triangle);

			MUVertex v = _mu_vertices[triangle.v0];
			v.tcount = 1;
			_mu_vertices.set(triangle.v0, v);

			v = _mu_vertices[triangle.v1];
			v.tcount = 1;
			_mu_vertices.set(triangle.v1, v);

			v = _mu_vertices[triangle.v2];
			v.tcount = 1;
			_mu_vertices.set(triangle.v2, v);
		}
	}

	_mu_triangles.resize(dst);

	dst = 0;
	for (int i = 0; i < _mu_vertices.size(); i++) {
		MUVertex vert = _mu_vertices[i];

		if (vert.tcount > 0) {
			vert.tstart = dst;
			_mu_vertices.set(i, vert);

			if (dst != i) {
				MUVertex dv = _mu_vertices[dst];
				dv.vertex = vert.vertex;
				_mu_vertices.set(dst, dv);

				if (_indices.size() > 0) _indices.set(dst, _indices[i]);
			}

			++dst;
		}
	}

	for (int i = 0; i < _mu_triangles.size(); i++) {
		MUTriangle triangle = _mu_triangles[i];
		triangle.v0 = _mu_vertices[triangle.v0].tstart;
		triangle.v1 = _mu_vertices[triangle.v1].tstart;
		triangle.v2 = _mu_vertices[triangle.v2].tstart;
		_mu_triangles.set(i, triangle);
	}

	//vertexCount = dst;
	if (_indices.size() > 0) _indices.resize(dst);
}

bool FastQuadraticMeshSimplifier::are_uvs_the_same(int channel, int indexA, int indexB) {
	Vector2 uva = _mu_vertices[indexA].vertex.uv;
	Vector2 uvb = _mu_vertices[indexB].vertex.uv;

	return Math::is_equal_approx(uva.x, uvb.x) && Math::is_equal_approx(uva.y, uvb.y);
}

/// Remove vertices and mark deleted triangles
void FastQuadraticMeshSimplifier::remove_vertex_pass(int startTrisCount, int targetTrisCount, double threshold, PoolVector<bool> *deleted0, PoolVector<bool> *deleted1, int *deletedTris) {
	Vector3 p;
	Vector3 barycentricCoord;
	for (int tid = 0; tid < _mu_triangles.size(); ++tid) {
		MUTriangle t = _mu_triangles[tid];

		if (t.dirty || t.deleted || t.err3 > threshold)
			continue;

		Vector3 errors = t.GetErrors();
		Vector3 attrib_indices = t.GetAttributeIndices();
		for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++) {
			if (errors[edgeIndex] > threshold)
				continue;

			int nextEdgeIndex = ((edgeIndex + 1) % 3);
			int i0 = t.get(edgeIndex);
			int i1 = t.get(nextEdgeIndex);

			MUVertex v0 = _mu_vertices[i0];
			MUVertex v1 = _mu_vertices[i1];

			// Border check
			if (v0.border_edge != v1.border_edge)
				continue;

			// Seam check
			else if (v0.uv_seam_edge != v1.uv_seam_edge)
				continue;

			// Foldover check
			else if (v0.uv_foldover_edge != v1.uv_foldover_edge)
				continue;

			// If borders should be preserved
			else if (_preserve_border_dges && v0.border_edge)
				continue;

			// If seams should be preserved
			else if (_preserve_uv_seam_edges && v0.uv_seam_edge)
				continue;

			// If foldovers should be preserved
			else if (_preserve_uv_foldover_edges && v0.uv_foldover_edge)
				continue;

			// Compute vertex to collapse to
			calculate_error(v0, v1, &p);
			deleted0->resize(v0.tcount); // normals temporarily
			deleted1->resize(v1.tcount); // normals temporarily

			// Don't remove if flipped
			if (flipped(p, i0, i1, v0, deleted0))
				continue;

			if (flipped(p, i1, i0, v1, deleted1))
				continue;

			// Calculate the barycentric coordinates within the triangle
			int nextNextEdgeIndex = ((edgeIndex + 2) % 3);
			int i2 = t.get(nextNextEdgeIndex);
			barycentricCoord = calculate_barycentric_coords(p, v0.vertex.vertex, v1.vertex.vertex, _mu_vertices[i2].vertex.vertex);

			// Not flipped, so remove edge
			v0.vertex.vertex = p;
			v0.q += v1.q;

			// Interpolate the vertex attributes
			int ia0 = attrib_indices[edgeIndex];
			int ia1 = attrib_indices[nextEdgeIndex];
			int ia2 = attrib_indices[nextNextEdgeIndex];
			interpolate_vertex_attributes(ia0, ia0, ia1, ia2, barycentricCoord);

			if (v0.uv_seam_edge) {
				ia0 = -1;
			}

			int tstart = _mu_refs.size();
			update_triangles(i0, ia0, v0, deleted0, deletedTris);
			update_triangles(i0, ia0, v1, deleted1, deletedTris);

			int tcount = _mu_refs.size() - tstart;
			if (tcount <= v0.tcount) {
				// save ram
				if (tcount > 0) {
					int dests = v0.tstart;
					for (int v = 0; v < tcount; ++v) {
						_mu_refs.set(v + tstart, _mu_refs[v + dests]);
					}
				}
			} else {
				// append
				//MUVertex v = _mu_vertices[i0];
				v0.tstart = tstart;
				//_mu_vertices.set(i0, v0);
			}

			//MUVertex v = _mu_vertices[i0];
			v0.tcount = tcount;

			_mu_vertices.set(i0, v0);
			break;
		}

		// Check if we are already done
		if ((startTrisCount - (*deletedTris)) <= targetTrisCount)
			break;
	}

	//return deletedTris;
}

double FastQuadraticMeshSimplifier::vertex_error(const SymmetricMatrix &q, const double x, const double y, const double z) const {
	return q.m0 * x * x + 2 * q.m1 * x * y + 2 * q.m2 * x * z + 2 * q.m3 * x + q.m4 * y * y + 2 * q.m5 * y * z + 2 * q.m6 * y + q.m7 * z * z + 2 * q.m8 * z + q.m9;
}

double FastQuadraticMeshSimplifier::calculate_error(const MUVertex &vert0, const MUVertex &vert1, Vector3 *result) {
	// compute interpolated vertex
	SymmetricMatrix q = (vert0.q + vert1.q);
	bool borderEdge = (vert0.border_edge & vert1.border_edge);
	double error = 0.0;
	double det = q.Determinant1();
	if (det != 0.0 && !borderEdge) {
		// q_delta is invertible
		result = new Vector3(
				-1.0 / det * q.Determinant2(), // vx = A41/det(q_delta)
				1.0 / det * q.Determinant3(), // vy = A42/det(q_delta)
				-1.0 / det * q.Determinant4()); // vz = A43/det(q_delta)
		error = vertex_error(q, result->x, result->y, result->z);
	} else {
		// det = 0 -> try to find best result
		Vector3 p1 = vert0.vertex.vertex;
		Vector3 p2 = vert1.vertex.vertex;
		Vector3 p3 = (p1 + p2) * 0.5f;
		double error1 = vertex_error(q, p1.x, p1.y, p1.z);
		double error2 = vertex_error(q, p2.x, p2.y, p2.z);
		double error3 = vertex_error(q, p3.x, p3.y, p3.z);

		error = FastQuadraticMeshSimplifier::min3(error1, error2, error3);
		if (error == error3) {
			result->x = p3.x;
			result->y = p3.y;
			result->z = p3.z;
		} else if (error == error2) {
			result->x = p2.x;
			result->y = p2.y;
			result->z = p2.z;
		} else if (error == error1) {
			result->x = p1.x;
			result->y = p1.y;
			result->z = p1.z;
		} else {
			result->x = p3.x;
			result->y = p3.y;
			result->z = p3.z;
		}
	}
	return error;
}

void FastQuadraticMeshSimplifier::update_triangles(int i0, int ia0, const MUVertex &v, PoolVector<bool> *deleted, int *deletedTriangles) {
	Vector3 p;
	int tcount = v.tcount;

	for (int k = 0; k < tcount; ++k) {
		MURef r = _mu_refs[v.tstart + k];
		int tid = r.tid;

		ERR_CONTINUE(r.tid >= _mu_triangles.size());

		MUTriangle t = _mu_triangles[tid];

		if (t.deleted)
			continue;

		if (deleted->get(k)) {
			t.deleted = true;

			_mu_triangles.set(tid, t);

			++(*deletedTriangles);
			continue;
		}

		t.set(r.tvertex, i0);
		if (ia0 != -1) {
			t.SetAttributeIndex(r.tvertex, ia0);
		}

		t.dirty = true;
		t.err0 = calculate_error(_mu_vertices[t.v0], _mu_vertices[t.v1], &p);
		t.err1 = calculate_error(_mu_vertices[t.v1], _mu_vertices[t.v2], &p);
		t.err2 = calculate_error(_mu_vertices[t.v2], _mu_vertices[t.v0], &p);
		t.err3 = FastQuadraticMeshSimplifier::min3(t.err0, t.err1, t.err2);

		_mu_triangles.set(tid, t);
		_mu_refs.push_back(r);
	}
}

bool FastQuadraticMeshSimplifier::flipped(const Vector3 &p, int i0, int i1, const MUVertex &v0, PoolVector<bool> *deleted) {
	int tcount = v0.tcount;

	for (int k = 0; k < tcount; k++) {
		MURef r = _mu_refs[v0.tstart + k];

		ERR_CONTINUE(r.tid >= _mu_triangles.size());

		MUTriangle t = _mu_triangles[r.tid];

		if (t.deleted)
			continue;

		int s = r.tvertex;
		int id1 = t.get((s + 1) % 3);
		int id2 = t.get((s + 2) % 3);
		if (id1 == i1 || id2 == i1) {
			deleted->set(k, true);
			continue;
		}

		Vector3 d1 = _mu_vertices[id1].vertex.vertex - p;
		d1.normalize();
		Vector3 d2 = _mu_vertices[id2].vertex.vertex - p;
		d2.normalize();
		double dot = d1.dot(d2);
		if (Math::abs(dot) > 0.999)
			return true;

		Vector3 n = d1.cross(d2);
		n.normalize();
		deleted->set(k, false);
		dot = n.dot(t.n);
		if (dot < 0.2)
			return true;
	}

	return false;
}

Vector3 FastQuadraticMeshSimplifier::calculate_barycentric_coords(Vector3 const &point, Vector3 const &a, Vector3 const &b, Vector3 const &c) {
	Vector3 v0 = (Vector3)(b - a), v1 = (Vector3)(c - a), v2 = (Vector3)(point - a);
	float d00 = v0.dot(v0);
	float d01 = v0.dot(v1);
	float d11 = v1.dot(v1);
	float d20 = v2.dot(v0);
	float d21 = v2.dot(v1);
	float denom = d00 * d11 - d01 * d01;
	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	float u = 1.0 - v - w;

	return Vector3(u, v, w);
}

void FastQuadraticMeshSimplifier::interpolate_vertex_attributes(int dst, int i0, int i1, int i2, const Vector3 &barycentricCoord) {
	MUVertex v0 = _mu_vertices[i0];
	MUVertex v1 = _mu_vertices[i1];
	MUVertex v2 = _mu_vertices[i2];
	MUVertex vdst = _mu_vertices[dst];

	if ((_format & VisualServer::ARRAY_FORMAT_NORMAL) != 0)
		vdst.vertex.normal = (v0.vertex.normal * barycentricCoord.x) + (v1.vertex.normal * barycentricCoord.y) + (v2.vertex.normal * barycentricCoord.z).normalized();

	if ((_format & VisualServer::ARRAY_FORMAT_TEX_UV) != 0)
		vdst.vertex.uv = (v0.vertex.uv * barycentricCoord.x) + (v1.vertex.uv * barycentricCoord.y) + (v2.vertex.uv * barycentricCoord.z);

	if ((_format & VisualServer::ARRAY_FORMAT_TEX_UV2) != 0)
		vdst.vertex.uv2 = (v0.vertex.uv2 * barycentricCoord.x) + (v1.vertex.uv2 * barycentricCoord.y) + (v2.vertex.uv2 * barycentricCoord.z);

	if ((_format & VisualServer::ARRAY_FORMAT_COLOR) != 0)
		vdst.vertex.color = (v0.vertex.color * barycentricCoord.x) + (v1.vertex.color * barycentricCoord.y) + (v2.vertex.color * barycentricCoord.z);

	_mu_vertices.set(dst, vdst);
}

void FastQuadraticMeshSimplifier::remove_doubles() {
	if (_vertices.size() == 0)
		return;

	//print_error("before " + String::num(_vertices.size()));

	for (int i = 0; i < _vertices.size(); ++i) {
		Vertex vert = _vertices[i];
		PoolVector<int> indices;

		for (int j = i + 1; j < _vertices.size(); ++j) {
			if (_vertices[j] == vert) {
				indices.push_back(j);
			}
		}

		for (int j = 0; j < indices.size(); ++j) {
			int index = indices[j];

			_vertices.remove(index);

			//make all indices that were bigger than the one we replaced one lower
			for (int k = 0; k < _indices.size(); ++k) {
				int indx = _indices[k];

				if (indx == index) {
					_indices.set(k, i);
				} else if (indx > index) {
					_indices.set(k, --indx);
				}
			}

			for (int k = j + 1; k < indices.size(); ++k) {
				int val = indices[k];

				if (val > index) {
					indices.set(k, --val);
				}
			}
		}
	}

	//print_error("after " + String::num(_vertices.size())+ " " + String::num(duration.count()));
}

//lot faster that normal remove_doubles, but false positives can happen curtesy of hash collisions
void FastQuadraticMeshSimplifier::remove_doubles_hashed() {
	if (_vertices.size() == 0)
		return;

	//print_error("before " + String::num(_vertices.size()));

	PoolVector<uint32_t> hashes;
	hashes.resize(_vertices.size());
	for (int i = 0; i < _vertices.size(); ++i) {
		hashes.set(i, VertexHasher::hash(_vertices[i]));
	}

	for (int i = 0; i < hashes.size(); ++i) {
		uint32_t hash = hashes[i];
		PoolVector<int> indices;

		for (int j = i + 1; j < hashes.size(); ++j) {
			if (hashes[j] == hash) {
				indices.push_back(j);
			}
		}

		for (int j = 0; j < indices.size(); ++j) {
			int index = indices[j];

			hashes.remove(index);
			_vertices.remove(index);

			//make all indices that were bigger than the one we replaced one lower
			for (int k = 0; k < _indices.size(); ++k) {
				int indx = _indices[k];

				if (indx == index) {
					_indices.set(k, i);
				} else if (indx > index) {
					_indices.set(k, --indx);
				}
			}

			for (int k = j + 1; k < indices.size(); ++k) {
				int val = indices[k];

				if (val > index) {
					indices.set(k, --val);
				}
			}
		}
	}

	//print_error("after " + String::num(_vertices.size()) + " " + String::num(duration.count()));
}

FastQuadraticMeshSimplifier::FastQuadraticMeshSimplifier() {
	_max_iteration_count = 100;
	_agressiveness = 7.0;
	_enable_smart_link = true;
	_preserve_border_dges = false;
	_preserve_uv_seam_edges = false;
	_preserve_uv_foldover_edges = false;
	_format = 0;
}

FastQuadraticMeshSimplifier::~FastQuadraticMeshSimplifier() {
}

void FastQuadraticMeshSimplifier::_bind_methods() {
	ClassDB::bind_method(D_METHOD("initialize", "arrays"), &FastQuadraticMeshSimplifier::initialize);
	ClassDB::bind_method(D_METHOD("get_arrays"), &FastQuadraticMeshSimplifier::get_arrays);
	ClassDB::bind_method(D_METHOD("simplify_mesh", "quality"), &FastQuadraticMeshSimplifier::simplify_mesh);
	ClassDB::bind_method(D_METHOD("simplify_mesh_lossless"), &FastQuadraticMeshSimplifier::simplify_mesh_lossless);

	ClassDB::bind_method(D_METHOD("get_max_iteration_count"), &FastQuadraticMeshSimplifier::get_max_iteration_count);
	ClassDB::bind_method(D_METHOD("set_max_iteration_count", "value"), &FastQuadraticMeshSimplifier::set_max_iteration_count);
	ADD_PROPERTY(PropertyInfo(Variant::INT, "max_iteration_count"), "set_max_iteration_count", "get_max_iteration_count");

	ClassDB::bind_method(D_METHOD("get_agressiveness"), &FastQuadraticMeshSimplifier::get_agressiveness);
	ClassDB::bind_method(D_METHOD("set_agressiveness", "value"), &FastQuadraticMeshSimplifier::set_agressiveness);
	ADD_PROPERTY(PropertyInfo(Variant::REAL, "agressiveness"), "set_agressiveness", "get_agressiveness");

	ClassDB::bind_method(D_METHOD("get_enable_smart_link"), &FastQuadraticMeshSimplifier::get_enable_smart_link);
	ClassDB::bind_method(D_METHOD("set_enable_smart_link", "value"), &FastQuadraticMeshSimplifier::set_enable_smart_link);
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "enable_smart_link"), "set_enable_smart_link", "get_enable_smart_link");

	ClassDB::bind_method(D_METHOD("get_preserve_border_dges"), &FastQuadraticMeshSimplifier::get_preserve_border_dges);
	ClassDB::bind_method(D_METHOD("set_preserve_border_dges", "value"), &FastQuadraticMeshSimplifier::set_preserve_border_dges);
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "preserve_border_dges"), "set_preserve_border_dges", "get_preserve_border_dges");

	ClassDB::bind_method(D_METHOD("get_preserve_uv_seam_edges"), &FastQuadraticMeshSimplifier::get_preserve_uv_seam_edges);
	ClassDB::bind_method(D_METHOD("set_preserve_uv_seam_edges", "value"), &FastQuadraticMeshSimplifier::set_preserve_uv_seam_edges);
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "preserve_uv_seam_edges"), "set_preserve_uv_seam_edges", "get_preserve_uv_seam_edges");

	ClassDB::bind_method(D_METHOD("get_preserve_uv_foldover_edges"), &FastQuadraticMeshSimplifier::get_preserve_uv_foldover_edges);
	ClassDB::bind_method(D_METHOD("set_preserve_uv_foldover_edges", "value"), &FastQuadraticMeshSimplifier::set_preserve_uv_foldover_edges);
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "preserve_uv_foldover_edges"), "set_preserve_uv_foldover_edges", "get_preserve_uv_foldover_edges");

	ClassDB::bind_method(D_METHOD("remove_doubles"), &FastQuadraticMeshSimplifier::remove_doubles);
	ClassDB::bind_method(D_METHOD("remove_doubles_hashed"), &FastQuadraticMeshSimplifier::remove_doubles_hashed);
}
