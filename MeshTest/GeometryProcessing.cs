using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using ShoNS.Array;
namespace mikity.GeometryProcessing
{
    class tuple
    {
        public face f;
        public int J;
        public tuple(face _f, int _J)
        {
            f = _f;
            J = _J;

        }
    }
    class face
    {
        public face(int _N, params int[] indices)
        {
            corner = indices;
            N = _N;
        }
        public int N;
        public int[] corner;
        public halfedge firsthalfedge;
    }

    class halfedge
    {
        public vertex P;
        public face owner;
        public halfedge pair, next, prev;
        //Warning
        //After boundary halfedges are appropriately setup, Naked halfedges should not exist.
        public bool isNaked
        {
            get
            {
                return pair == null ? true : false;
            }
        }
        public bool isBoundary
        {
            get
            {
                return owner == null ? true : false;
            }
        }
        public halfedge(vertex _P)
        {
            P = _P;
            if (_P.hf_begin == null) _P.hf_begin = this;
        }
    }
    class vertex
    {
        public int N;
        public List<halfedge> star = new List<halfedge>();
        public List<halfedge> onering = new List<halfedge>();
        public halfedge hf_begin;
        //public halfedge hf_end;
        /*public bool isNaked
        {
            get
            {
                return hf_begin == hf_end ? false : true;
            }
        }*/
        public vertex(int _N)
        {
            N = _N;
        }
        public bool isInner
        {
            get
            {
                return onering[0] == onering[onering.Count - 1].next;
            }
        }
        public bool isBoundary
        {
            get
            {
                return onering[0] != onering[onering.Count - 1].next;
            }
        }
    }
    class MeshStructure
    {
        //to get boundary chain
        //boundaryStart->hf->next->P->hf->next->P->....
        //
        //public vertex boundaryStart;
        public halfedge[] boundaryStart;
        public List<vertex> vertices = new List<vertex>();
        public List<face> faces = new List<face>();
        public List<halfedge> halfedges = new List<halfedge>();
        public List<vertex> innerVertices = new List<vertex>();
        public List<vertex> outerVertices = new List<vertex>();  
        //private halfedge[,] __halfedgeTable;
        //private List<face>[,] _faceTable;
        List<halfedge>[] __halfedgeTable;
        List<tuple>[] _faceTable;
        private orient[] __orientation;
        public List<halfedge> edges()
        {
            var res = new List<halfedge>();
            foreach (var e in halfedges)
            {
                if (e.P.N < e.next.P.N)
                    res.Add(e);
            }
            return res;
        }
        public int nVertices
        {
            get
            {
                return vertices.Count;
            }
        }
        public int nFaces
        {
            get
            {
                return faces.Count;
            }
        }
        private MeshStructure()
        {
            this.Clear();
        }
        private enum orient
        {
            unknown, positive, negative
        }
        void getfaces(int I, int J, List<face> faces)
        {
            if (faces == null) faces = new List<face>();
            faces.Clear();
            var ff = _faceTable[I];
            foreach (var f in ff)
            {
                if (f.J==J) faces.Add(f.f);
            }
        }
        void gethalfedges(int I, int J, List<halfedge> lhalfedges)
        {
            if (lhalfedges == null) lhalfedges  = new List<halfedge>();
            lhalfedges.Clear();
            var ff = __halfedgeTable[I];
            foreach (var f in ff)
            {
                if (f.next.P.N == J) lhalfedges.Add(f);
            }
        }

        private void Construct(Mesh val)
        {
            int _nVertices = val.Vertices.Count;
            int _nFaces = val.Faces.Count;

            __orientation = new orient[_nFaces];
            //_faceTable = new List<face>[_nVertices, _nVertices];
            //__halfedgeTable = new halfedge[_nVertices, _nVertices];
            
            _faceTable = new List<tuple>[_nVertices];
            __halfedgeTable = new List<halfedge>[_nVertices];
            for (int i = 0; i < _nVertices; i++)
            {
                _faceTable[i] = new List<tuple>();
                __halfedgeTable[i] = new List<halfedge>();
            }
            for (int i = 0; i < __orientation.Count(); i++)
            {
                __orientation[i] = orient.unknown;
            }

            for (int i = 0; i < _nVertices; i++)
            {
                var _v = new vertex(i);
                vertices.Add(_v);
            }

            for (int i = 0; i < _nFaces; i++)
            {
                var f = val.Faces[i];
                var _f = f.IsQuad ? new face(i, f.A, f.B, f.C, f.D) : new face(i, f.A, f.B, f.C);
                faces.Add(_f);
                faceTableAdd(_f);
            }
            //Recursive
            halfEdgeAdd(faces[0]);
            List<halfedge> lhalfedges = new List<halfedge>();
            //find pairs
            foreach (var h in halfedges)
            {
                int i = h.P.N;
                int j = h.next.P.N;
                gethalfedges(i, j, lhalfedges);
                if (lhalfedges.Count!=0) throw new ArgumentOutOfRangeException(";)");
                //__halfedgeTable[i, j] = h;
                __halfedgeTable[i].Add(h);
            }
            foreach (var h in halfedges)
            {
                int i = h.P.N;
                int j = h.next.P.N;
                gethalfedges(j, i, lhalfedges);
                //if boundary edge...
                if (lhalfedges.Count==0)
                {
                    h.pair = null;
                }
                else
                {
                    h.pair = lhalfedges[0];
                }
            }
            //post process to find boundary vertices

            //align the first half edge at boundary
            foreach (var v in vertices)
            {
                var h = v.hf_begin;
                //v.hf_end = h;
                do
                {
                    if (h.prev.isNaked)
                    {
                        //v.hf_end = h.prev;
                        while (!h.isNaked)
                        {
                            h = h.pair.next;
                        }
                        v.hf_begin = h;
                        break;
                    }
                    h = h.prev.pair;
                } while (h != v.hf_begin);
            }
            List<List<halfedge>> boundary_complements = new List<List<halfedge>>();
            int nBoundary = 0;
            foreach (var v in vertices)
            {
                var h = v.hf_begin;
                if (h.isNaked)//first naked halfedge found
                {
                    bool flag = true;
                    for (int i = 0; i < nBoundary; i++)
                    {
                        if(boundary_complements[i].Contains(h))
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag)
                    {
                        boundary_complements.Add(new List<halfedge>());
                        do
                        {
                            boundary_complements[nBoundary].Add(h);
                            h = h.next.P.hf_begin;
                        } while (h != v.hf_begin);
                        nBoundary++;
                    }
                    //break;
                }
            }
            //insert boundary halfedges
            boundaryStart = new halfedge[nBoundary];
            foreach (var boundary_complement in boundary_complements)
            {
                List<halfedge> boundary = new List<halfedge>();
                for (int i = 0; i < boundary_complement.Count; i++)
                {
                    boundary.Add(new halfedge(boundary_complement[i].next.P));
                    boundary[i].pair = boundary_complement[i];
                    boundary_complement[i].pair = boundary[i];
                }
                halfedges.AddRange(boundary);
                boundaryStart[boundary_complements.IndexOf(boundary_complement)] = boundary[0];
                for (int i = 0; i < boundary.Count; i++)
                {
                    boundary[i].owner = null;
                    if (i != 0)
                    {
                        boundary[i].next = boundary_complement[i - 1].pair;
                    }
                    else
                    {
                        boundary[i].next = boundary_complement[boundary_complement.Count - 1].pair;
                    }
                    if (i != boundary.Count - 1)
                    {
                        boundary[i].prev = boundary_complement[i + 1].pair;
                    }
                    else
                    {
                        boundary[i].prev = boundary_complement[0].pair;
                    }
                }
            }
            //check if any naked halfedge survives
            foreach (var e in halfedges)
            {
                if (e.isNaked) System.Windows.Forms.MessageBox.Show("error");
            }
            //post process to create stars
            foreach (var v in vertices)
            {
                var h = v.hf_begin;
                v.star.Clear();
                do
                {
                    v.star.Add(h);
                    if (h.isBoundary) break;
                    h = h.prev.pair;
                } while (h != v.hf_begin);
            }
            //post process to create onering
            foreach (var v in vertices)
            {
                var h = v.hf_begin;
                v.onering.Clear();
                do
                {
                    do
                    {
                        h = h.next;
                        v.onering.Add(h);
                    } while (h.next.next.P != v);
                    if (h.next.pair.isBoundary) break;
                    h = h.next.pair;
                } while (h != v.hf_begin);
            }
            //post process to split the vertices into inner and outer.
            innerVertices.Clear();
            outerVertices.Clear();
            foreach (var v in vertices)
            {
                if (v.hf_begin.pair.isBoundary) outerVertices.Add(v); else innerVertices.Add(v);
            }
        }

                
        private void halfEdgeAdd(face f)
        {
            var _o = orient.unknown;
            List<face> faces=new List<face>();
            for (int i = 0; i < f.corner.Count(); i++)
            {
                int I = f.corner[i];
                int J = (i == f.corner.Count() - 1) ? f.corner[0] : f.corner[i + 1];
                getfaces(I, J, faces);                
                if (faces.Count == 2)
                {
                    if (faces[0] == f)
                    {
                        if (__orientation[faces[1].N] != orient.unknown)
                        {
                            _o = __orientation[faces[1].N] == orient.positive ? orient.negative : orient.positive;
                        }
                    }
                    if (faces[1] == f)
                    {
                        if (__orientation[faces[0].N] != orient.unknown)
                        {
                            _o = __orientation[faces[0].N] == orient.positive ? orient.negative : orient.positive;
                        }
                    }
                }
                else
                {
                    getfaces(J, I, faces);
                    if (faces.Count != 0)
                    {
                        if (__orientation[faces[0].N] != orient.unknown)
                        {
                            _o = __orientation[faces[0].N];
                        }
                    }
                }
            }
            __orientation[f.N] = _o == orient.unknown ? orient.positive : _o;
            //register a halfedge
            if (f.corner.Count() == 3 && __orientation[f.N] == orient.positive)
            {
                var he1 = new halfedge(vertices[f.corner[0]]);
                var he2 = new halfedge(vertices[f.corner[1]]);
                var he3 = new halfedge(vertices[f.corner[2]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                he1.prev = he3; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he1; he3.owner = f;
                f.firsthalfedge = he1;
            }

            if (f.corner.Count() == 3 && __orientation[f.N] == orient.negative)
            {
                var he1 = new halfedge(vertices[f.corner[2]]);
                var he2 = new halfedge(vertices[f.corner[1]]);
                var he3 = new halfedge(vertices[f.corner[0]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                he1.prev = he3; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he1; he3.owner = f;
                f.firsthalfedge = he1;
            }

            if (f.corner.Count() == 4 && __orientation[f.N] == orient.positive)
            {
                var he1 = new halfedge(vertices[f.corner[0]]);
                var he2 = new halfedge(vertices[f.corner[1]]);
                var he3 = new halfedge(vertices[f.corner[2]]);
                var he4 = new halfedge(vertices[f.corner[3]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                halfedges.Add(he4);
                he1.prev = he4; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he4; he3.owner = f;
                he4.prev = he3; he4.next = he1; he4.owner = f;
                f.firsthalfedge = he1;
            }

            if (f.corner.Count() == 4 && __orientation[f.N] == orient.negative)
            {
                var he1 = new halfedge(vertices[f.corner[3]]);
                var he2 = new halfedge(vertices[f.corner[2]]);
                var he3 = new halfedge(vertices[f.corner[1]]);
                var he4 = new halfedge(vertices[f.corner[0]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                halfedges.Add(he4);
                he1.prev = he4; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he4; he3.owner = f;
                he4.prev = he3; he4.next = he1; he4.owner = f;
                f.firsthalfedge = he1;
            }

            //list up neighbors that are not oriented
            for (int i = 0; i < f.corner.Count(); i++)
            {
                int I = f.corner[i];
                int J = (i == f.corner.Count() - 1) ? f.corner[0] : f.corner[i + 1];
                getfaces(I, J, faces);
                if (faces.Count == 2)
                {
                    if (faces[0] == f)
                    {
                        if (__orientation[faces[1].N] == orient.unknown)
                        {
                            halfEdgeAdd(faces[1]);
                        }
                    }
                    if (faces[1] == f)
                    {
                        if (__orientation[faces[0].N] == orient.unknown)
                        {
                            halfEdgeAdd(faces[0]);
                        }
                    }
                }
                else
                {
                    getfaces(J, I, faces);
                    if (faces.Count != 0)
                    {
                        if (__orientation[faces[0].N] == orient.unknown)
                        {
                            halfEdgeAdd(faces[0]);
                        }
                    }
                }
            }
        }
        private void faceTableAdd(int i,int j,face f)
        {
            _faceTable[i].Add(new tuple(f,j));
        }
        private void faceTableAdd(face f)
        {
            for (int i = 0; i < f.corner.Count(); i++)
            {
                int I = f.corner[i];
                int J = (i == f.corner.Count() - 1) ? f.corner[0] : f.corner[i + 1];
                faceTableAdd(I, J,f);
            }
        }
        public static MeshStructure CreateFrom(Mesh val)
        {
            var ret = new MeshStructure();
            ret.Construct(val);
            return ret;
        }
        public void Clear()
        {
            vertices.Clear();
            faces.Clear();
            halfedges.Clear();
            innerVertices.Clear();
            outerVertices.Clear();
            boundaryStart = null;
        }
    }
}
