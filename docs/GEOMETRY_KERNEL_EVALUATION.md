# Geometry Kernel Evaluation

## 1. Executive Summary

**Recommendation:** We should use **Gmsh's built-in OpenCASCADE (OCC) kernel** for all procedural geometry generation.

**Reasoning:** This approach offers the most direct and robust path to our end-to-end simulation goal. It introduces **zero new dependencies** into our already-working solver environment, radically simplifies the project's architecture, and provides an API perfectly suited for a headless, automated pipeline. While other tools are more powerful as general-purpose CAD applications, Gmsh is the optimal choice for this specific, script-driven simulation task.

## 2. The Core Problem: Headless, Containerized Automation

Our goal is not to create a desktop CAD application, but an automated, containerized pipeline that runs without human intervention. This context is critical. The primary challenge is not "Which tool has the most features?" but "Which tool is best suited for a lightweight, scriptable, and reliable service that generates geometry and passes it to a solver?"

## 3. Evaluation of Alternatives

### Option A: FreeCAD (The Initial Idea)

FreeCAD is a powerful, full-featured parametric 3D CAD modeler with an extensive Python API.

-   **Dependency Overhead: Very High.** FreeCAD is a large desktop application. To run it headlessly in a container requires installing numerous dependencies, including GUI libraries (even if unused), and often requires complex entrypoint scripts involving Xvfb (a virtual X-server). This would result in a very large and complex Docker image for the `horn-geometry` service.
-   **Architectural Simplicity: Poor.** This approach necessitates a separate, heavy `horn-geometry` service. The pipeline would involve one large container building a STEP file, saving it to a shared volume, and a second large container picking it up for meshing and solving. This increases complexity and potential points of failure.
-   **API Suitability: Good, but...** The FreeCAD Python API is powerful but is fundamentally designed to script the application, not as a standalone library.
-   **Meshing Integration: Poor.** The integration is entirely indirect, via the export/import of a STEP file.

### Option B: Gmsh's Internal CAD Kernel (The Recommended Approach)

Gmsh is not just a mesher; it includes its own full-featured geometry kernel based on OpenCASCADE, the same technology that underlies FreeCAD.

-   **Dependency Overhead: None.** This is the decisive advantage. The `dolfinx/dolfinx` base image **already includes a full installation of Gmsh with its Python API.** We do not need to install anything new. The dependency is already met.
-   **Architectural Simplicity: Excellent.** We do not need a separate `horn-geometry` service. The geometry generation can be a simple Python function within the existing `horn-solver` package or a very lightweight `horn-geometry` package that has no heavy dependencies. The entire "geometry -> mesh" pipeline can happen within a single script and even in-memory, without ever writing a STEP file to disk.
-   **API Suitability: Excellent.** The Gmsh Python API is explicitly designed for scripting and headless automation. It is lightweight and robust.
-   **Meshing Integration: Perfect.** This is the key workflow benefit. We can build the geometry and mesh it in the same script, using the same tool. We can directly reference geometric entities during meshing, which is far more robust than relying on the blind import of a STEP file.

### Option C: Other Python-Native Libraries (e.g., CadQuery)

Libraries like CadQuery are designed from the ground up to be "CAD-as-code" and are excellent alternatives.

-   **Dependency Overhead: Low.** These are standard Python packages that can be installed with `pip`. While they have dependencies on OCC bindings (`pyocct`), they are far lighter than installing the entire FreeCAD application.
-   **Architectural Simplicity: Good.** Similar to the Gmsh approach, this would allow for a lightweight `horn-geometry` service.
-   **API Suitability: Excellent.** CadQuery, in particular, has a very clean, modern, and fluent API for programmatic CAD.
-   **Meshing Integration: Poor.** Like FreeCAD, the integration with our solver would be indirect, requiring the export of a STEP file.

## 4. Final Conclusion

| Criteria                   | FreeCAD      | Gmsh (Internal Kernel) | CadQuery       |
| -------------------------- | ------------ | ---------------------- | -------------- |
| **Dependency Overhead**    | Very High    | **None**               | Low            |
| **Architectural Simplicity** | Poor         | **Excellent**          | Good           |
| **API Suitability**        | Good         | **Excellent**          | Excellent      |
| **Meshing Integration**    | Poor         | **Perfect**            | Poor           |
| **Path to E2E MVP**        | Slow/Complex | **Fastest/Simplest**   | Medium         |

While CadQuery is an excellent modern library, the fact that **Gmsh is already a required, installed dependency for our meshing step** makes it the overwhelmingly logical choice. Using its internal kernel provides the fastest, simplest, and most robust path forward to achieving our goal of a minimal end-to-end solution. It eliminates an entire layer of architectural complexity and potential dependency issues. 