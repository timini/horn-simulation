import os
import psycopg2

def get_driver_parameters(driver_id: str) -> dict:
    """
    Retrieves Thiele/Small parameters for a given driver from the database.

    Args:
        driver_id: The unique identifier for the driver.

    Returns:
        A dictionary containing the driver's T/S parameters.
    
    Raises:
        ValueError: If the driver_id is not found in the database.
    """
    try:
        conn = psycopg2.connect(
            dbname=os.getenv("POSTGRES_DB", "horn_db"),
            user=os.getenv("POSTGRES_USER", "horn_user"),
            password=os.getenv("POSTGRES_PASSWORD", "horn_password"),
            host=os.getenv("POSTGRES_HOST", "localhost"),
            port=os.getenv("POSTGRES_PORT", "5432"),
        )
    except psycopg2.OperationalError as e:
        raise ConnectionError("Could not connect to the database. Is it running?") from e

    with conn.cursor() as cur:
        cur.execute("SELECT * FROM drivers WHERE driver_id = %s", (driver_id,))
        
        row = cur.fetchone()
        if row is None:
            conn.close()
            raise ValueError(f"Driver with ID '{driver_id}' not found in the database.")

        # Fetch column names from the cursor description to build a dictionary
        colnames = [desc[0] for desc in cur.description]
        driver_params = dict(zip(colnames, row))
    
    conn.close()
    
    return driver_params 