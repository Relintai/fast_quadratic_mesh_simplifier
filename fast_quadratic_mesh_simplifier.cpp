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

void FastQuadraticMeshSimplifier::initialize(const Array &arrays) {
	ERR_FAIL_COND(arrays.size() != ArrayMesh::ARRAY_MAX);

	PoolVector<Vector3> vertices = arrays.get(ArrayMesh::ARRAY_VERTEX);
	_vertices.resize(vertices.size());
	for (int i = 0; i < vertices.size(); ++i) {
		_vertices.set(i, vertices[i]);
	}

	PoolVector<Vector3> normals = arrays.get(ArrayMesh::ARRAY_NORMAL);
	_normals.resize(normals.size());
	for (int i = 0; i < normals.size(); ++i) {
		_normals.set(i, normals[i]);
	}

	PoolVector<Color> colors = arrays.get(ArrayMesh::ARRAY_COLOR);
	_colors.resize(colors.size());
	for (int i = 0; i < colors.size(); ++i) {
		_colors.set(i, colors[i]);
	}

	PoolVector<Vector2> uvs = arrays.get(ArrayMesh::ARRAY_TEX_UV);
	_uvs.resize(uvs.size());
	for (int i = 0; i < uvs.size(); ++i) {
		_uvs.set(i, uvs[i]);
	}

	PoolVector<Vector2> uv2s = arrays.get(ArrayMesh::ARRAY_TEX_UV2);
	_uv2s.resize(uv2s.size());
	for (int i = 0; i < uv2s.size(); ++i) {
		_uv2s.set(i, uv2s[i]);
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

void FastQuadraticMeshSimplifier::refresh_vertices() {
	_vertices.resize(_mu_vertices.size());
	for (int i = 0; i < _mu_vertices.size(); ++i) {
		MUVertex vert = _mu_vertices[i];

		_vertices.set(i, Vector3(vert.p));
	}
}

Array FastQuadraticMeshSimplifier::get_arrays() {
	Array arr;

	arr.resize(ArrayMesh::ARRAY_MAX);

	arr.set(ArrayMesh::ARRAY_VERTEX, _vertices);
	arr.set(ArrayMesh::ARRAY_NORMAL, _normals);
	arr.set(ArrayMesh::ARRAY_COLOR, _colors);
	arr.set(ArrayMesh::ARRAY_TEX_UV, _uvs);
	arr.set(ArrayMesh::ARRAY_TEX_UV2, _uv2s);
	arr.set(ArrayMesh::ARRAY_INDEX, _indices);

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

//Mesh Simplification
//Ported from https://github.com/Whinarn/UnityFastQuadraticMeshSimplifier
//Original license: MIT License Copyright (c) 2017 Mattias Edlund
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
			if (!_mu_triangles[i].deleted) {
				if (dst != i) {
					_mu_triangles[dst] = _mu_triangles[i];
				}
				dst++;
			}
		}
		_mu_triangles.resize(dst);
	}

	update_references();

	// Identify boundary : vertices[].border=0,1
	if (iteration == 0) {
		PoolVector<int> vcount;
		vcount.resize(8);
		PoolVector<int> vids;
		vids.resize(8);

		int vsize = 0;
		for (int i = 0; i < _mu_vertices.size(); i++) {
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
			int tstart = _mu_vertices[i].tstart;
			int tcount = _mu_vertices[i].tcount;
			vcount.resize(0);
			vids.resize(0);
			vsize = 0;

			for (int j = 0; j < tcount; j++) {
				int tid = _mu_refs[tstart + j].tid;
				MUTriangle t = _mu_triangles[tid];

				for (int k = 0; k < 3; k++) {
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

			for (int j = 0; j < vsize; j++) {
				if (vcount[j] == 1) {
					id = vids[j];

					MUVertex v = _mu_vertices[id];

					v.border_edge = true;
					_mu_vertices.set(id, v);

					++borderVertexCount;

					if (_enable_smart_link) {
						if (v.p.x < borderMinX) {
							borderMinX = v.p.x;
						}
						if (v.p.x > borderMaxX) {
							borderMaxX = v.p.x;
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
			for (int i = 0; i < _mu_vertices.size(); i++) {
				if (_mu_vertices[i].border_edge) {
					int vertexHash = (int)(((((_mu_vertices[i].p.x - borderMinX) / borderAreaWidth) * 2.0) - 1.0) * std::numeric_limits<int>::max());
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
			for (int i = 0; i < borderIndexCount; i++) {
				int myIndex = borderVertices[i].index;
				if (myIndex == -1)
					continue;

				Vector3 myPoint = _mu_vertices[myIndex].p;
				for (int j = i + 1; j < borderIndexCount; j++) {
					int otherIndex = borderVertices[j].index;
					if (otherIndex == -1)
						continue;
					else if ((borderVertices[j].hash - borderVertices[i].hash) > hashMaxDistance) // There is no point to continue beyond this point
						break;

					Vector3 otherPoint = _mu_vertices[otherIndex].p;
					double sqrX = ((myPoint.x - otherPoint.x) * (myPoint.x - otherPoint.x));
					double sqrY = ((myPoint.y - otherPoint.y) * (myPoint.y - otherPoint.y));
					double sqrZ = ((myPoint.z - otherPoint.z) * (myPoint.z - otherPoint.z));
					double sqrMagnitude = sqrX + sqrY + sqrZ;

					if (sqrMagnitude <= _vertex_link_distance_sqr) {
						borderVertices.get(j).set_index(-1); // NOTE: This makes sure that the "other" vertex is not processed again

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
			_mu_vertices[i].q.reset();
		}

		int v0, v1, v2;
		Vector3 n, p0, p1, p2, p10, p20, dummy;
		SymmetricMatrix sm;
		for (int i = 0; i < _mu_triangles.size(); ++i) {
			v0 = _mu_triangles[i].v0;
			v1 = _mu_triangles[i].v1;
			v2 = _mu_triangles[i].v2;

			p0 = _mu_vertices[v0].p;
			p1 = _mu_vertices[v1].p;
			p2 = _mu_vertices[v2].p;
			p10 = p1 - p0;
			p20 = p2 - p0;

			n = p10.cross(p20);

			n.normalize();
			_mu_triangles[i].n = n;

			sm.from_plane(n.x, n.y, n.z, -n.dot(p0));
			_mu_vertices[v0].q += sm;
			_mu_vertices[v1].q += sm;
			_mu_vertices[v2].q += sm;
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
		v.tcount++;
		_mu_vertices.set(t.v0, v);

		v = _mu_vertices[t.v1];
		v.tcount++;
		_mu_vertices.set(t.v1, v);

		v = _mu_vertices[t.v2];
		v.tcount++;
		_mu_vertices.set(t.v2, v);
	}

	int tstart = 0;
	for (int i = 0; i < _mu_vertices.size(); i++) {
		MUVertex v = _mu_vertices[i];

		v.tstart = tstart;
		tstart += v.tcount;
		v.tcount = 0;
		_mu_vertices.set(i, v);
	}

	// Write References
	_mu_refs.resize(tstart);
	for (int i = 0; i < _mu_triangles.size(); i++) {
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

	for (int i = 0; i < _mu_vertices.size(); i++) {
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
				d.p = s.p;
				_mu_vertices.set(iDest, d);

				triangle.v0 = triangle.va0;
			}

			if (triangle.va1 != triangle.v1) {
				int iDest = triangle.va1;
				int iSrc = triangle.v1;

				MUVertex d = _mu_vertices[iDest];
				MUVertex s = _mu_vertices[iSrc];
				d.p = s.p;
				_mu_vertices.set(iDest, d);

				triangle.v1 = triangle.va1;
			}

			if (triangle.va2 != triangle.v2) {
				int iDest = triangle.va2;
				int iSrc = triangle.v2;

				MUVertex d = _mu_vertices[iDest];
				MUVertex s = _mu_vertices[iSrc];
				d.p = s.p;
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
			_mu_vertices[i] = vert;

			if (dst != i) {
				_mu_vertices[dst].p = vert.p;

				if (_normals.size() > 0) _normals.set(dst, _normals[i]);

				if (_colors.size() > 0) _colors.set(dst, _colors[i]);
				if (_uvs.size() > 0) _uvs.set(dst, _uvs[i]);
				if (_uv2s.size() > 0) _uv2s.set(dst, _uv2s[i]);
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
		_mu_triangles[i] = triangle;
	}

	//vertexCount = dst;
	_vertices.resize(dst);
	if (_normals.size() > 0) _normals.resize(dst);
	if (_colors.size() > 0) _colors.resize(dst);
	if (_uvs.size() > 0) _uvs.resize(dst);
	if (_uv2s.size() > 0) _uv2s.resize(dst);
	if (_indices.size() > 0) _indices.resize(dst);
}

bool FastQuadraticMeshSimplifier::are_uvs_the_same(int channel, int indexA, int indexB) {
	if (_uv2s.size() > 0) {
		//Vector2 vertUV = _uv2s[channel];

		Vector2 uvA = _uv2s[indexA];
		Vector2 uvB = _uv2s[indexB];
		return uvA == uvB;
	}

	return false;
}

/// Remove vertices and mark deleted triangles
void FastQuadraticMeshSimplifier::remove_vertex_pass(int startTrisCount, int targetTrisCount, double threshold, PoolVector<bool> *deleted0, PoolVector<bool> *deleted1, int *deletedTris) {
	Vector3 p;
	Vector3 barycentricCoord;
	for (int tid = 0; tid < _mu_triangles.size(); tid++) {
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
			barycentricCoord = calculate_barycentric_coords(p, v0.p, v1.p, _mu_vertices[i2].p);

			// Not flipped, so remove edge
			v0.p = p;
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

double FastQuadraticMeshSimplifier::vertex_error(SymmetricMatrix q, double x, double y, double z) {
	return q.m0 * x * x + 2 * q.m1 * x * y + 2 * q.m2 * x * z + 2 * q.m3 * x + q.m4 * y * y + 2 * q.m5 * y * z + 2 * q.m6 * y + q.m7 * z * z + 2 * q.m8 * z + q.m9;
}

double FastQuadraticMeshSimplifier::calculate_error(MUVertex vert0, MUVertex vert1, Vector3 *result) {
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
		Vector3 p1 = vert0.p;
		Vector3 p2 = vert1.p;
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

	for (int k = 0; k < tcount; k++) {
		MURef r = _mu_refs[v.tstart + k];
		int tid = r.tid;
		MUTriangle t = _mu_triangles[tid];
		if (t.deleted)
			continue;

		if (deleted->get(k)) {
			MUTriangle t2 = _mu_triangles[tid];

			t2.deleted = true;

			_mu_triangles.set(tid, t2);

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

		_mu_triangles[tid] = t;
		_mu_refs.push_back(r);
	}

	//return deletedTriangles;
}

bool FastQuadraticMeshSimplifier::flipped(const Vector3 &p, int i0, int i1, const MUVertex &v0, PoolVector<bool> *deleted) {
	int tcount = v0.tcount;

	for (int k = 0; k < tcount; k++) {
		MURef r = _mu_refs[v0.tstart + k];
		MUTriangle t = _mu_triangles[r.tid];

		if (_mu_triangles[r.tid].deleted)
			continue;

		int s = r.tvertex;
		int id1 = t.get((s + 1) % 3);
		int id2 = t.get((s + 2) % 3);
		if (id1 == i1 || id2 == i1) {
			deleted->set(k, true);
			continue;
		}

		Vector3 d1 = _mu_vertices[id1].p - p;
		d1.normalize();
		Vector3 d2 = _mu_vertices[id2].p - p;
		d2.normalize();
		double dot = d1.dot(d2);
		if (Math::abs(dot) > 0.999)
			return true;

		Vector3 n = d1.cross(d2);
		n.normalize();
		deleted->set(k, false);
		dot = n.dot(_mu_triangles[r.tid].n);
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

void FastQuadraticMeshSimplifier::interpolate_vertex_attributes(int dst, int i0, int i1, int i2, Vector3 &barycentricCoord) {
	if (_normals.size() > 0) {
		_normals[dst] = (_normals[i0] * barycentricCoord.x) + (_normals[i1] * barycentricCoord.y) + (_normals[i2] * barycentricCoord.z).normalized();
	}

	if (_uvs.size() > 0) {
		_uvs[dst] = (_uvs[i0] * barycentricCoord.x) + (_uvs[i1] * barycentricCoord.y) + (_uvs[i2] * barycentricCoord.z);
	}

	if (_uv2s.size() > 0) {
		_uv2s[dst] = (_uv2s[i0] * barycentricCoord.x) + (_uv2s[i1] * barycentricCoord.y) + (_uv2s[i2] * barycentricCoord.z);
	}

	if (_colors.size() > 0) {
		_colors[dst] = (_colors[i0] * barycentricCoord.x) + (_colors[i1] * barycentricCoord.y) + (_colors[i2] * barycentricCoord.z);
	}
}

FastQuadraticMeshSimplifier::FastQuadraticMeshSimplifier() {
	_max_iteration_count = 100;
	_agressiveness = 7.0;
	_enable_smart_link = true;
	_preserve_border_dges = false;
	_preserve_uv_seam_edges = false;
	_preserve_uv_foldover_edges = false;
}

void FastQuadraticMeshSimplifier::_bind_methods() {
	ClassDB::bind_method(D_METHOD("initialize", "arrays"), &FastQuadraticMeshSimplifier::initialize);
	ClassDB::bind_method(D_METHOD("get_arrays"), &FastQuadraticMeshSimplifier::get_arrays);
	ClassDB::bind_method(D_METHOD("simplify_mesh", "quality"), &FastQuadraticMeshSimplifier::simplify_mesh);
	ClassDB::bind_method(D_METHOD("simplify_mesh_lossless"), &FastQuadraticMeshSimplifier::simplify_mesh_lossless);

	//	ClassDB::bind_method(D_METHOD("get_body_path"), &FastQuadraticMeshSimplifier::get_body_path);
	//	ClassDB::bind_method(D_METHOD("set_body_path", "value"), &FastQuadraticMeshSimplifier::set_body_path);
	//	ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "body_path"), "set_body_path", "get_body_path");
}