# 01: Project Setup and Database Integration

**Status:** Completed ✅

This task covered the initial setup of the project structure and the integration of a containerized PostgreSQL database.

## Key Accomplishments

- **Project Skeleton**:
    - Created a modular directory structure for source code (`src/horn`) and tests (`tests`).
    - Set up a `pyproject.toml` file using Poetry for dependency management.
    - Added a `README.md` with project information and a `.gitignore` file.
- **Test-Driven Development (TDD)**:
    - Established a TDD workflow with `pytest`.
    - Wrote an initial failing test for the (now-defunct) PDF parsing module and made it pass.
- **Database Setup**:
    - Created a `docker-compose.yml` file to run a PostgreSQL database in a Docker container.
    - Wrote an `init.sql` script to create the `drivers` table and seed it with a sample driver.
    - Implemented a `driver_db.py` module to handle database connections and fetch driver data.
- **Database Testing**:
    - Wrote tests for the `driver_db` module to ensure data could be retrieved successfully and errors were handled correctly.
    - Configured `pytest` to mark and run database-dependent tests. 