# dolfinx supplies native gmsh, including on Apple Silicon / Linux arm64.
FROM dolfinx/dolfinx:v0.8.0 AS production
WORKDIR /app
COPY ./packages/horn-core /app/packages/horn-core
COPY ./packages/horn-geometry /app/packages/horn-geometry
# Do not replace the system gmsh with an unavailable arm64 PyPI wheel.
RUN pip install --no-deps /app/packages/horn-core /app/packages/horn-geometry \
    && python3 -c "import gmsh, numpy; import horn_geometry.generator"

FROM production AS test
RUN pip install pytest
COPY ./packages/horn-geometry/tests /app/packages/horn-geometry/tests
