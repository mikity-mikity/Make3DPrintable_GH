using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mikity.ghComponents
{
    public partial class MeshTest : Grasshopper.Kernel.GH_Component
    {
        public override void BakeGeometry(Rhino.RhinoDoc doc, Rhino.DocObjects.ObjectAttributes att, List<Guid> obj_ids)
        {
            if (newMesh != null)
            {
                Rhino.DocObjects.ObjectAttributes a2 = att.Duplicate();
                a2.LayerIndex = 2;
                Guid id = doc.Objects.AddMesh(newMesh, a2);
                obj_ids.Add(id);
            }
        }
    }
}
