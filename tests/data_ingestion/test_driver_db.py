import pytest
from horn.data_ingestion import driver_db

# Mark all tests in this file as 'db'
pytestmark = pytest.mark.db

def test_get_driver_parameters_success():
    """
    Test retrieving an existing driver from the database.
    This test requires the Docker container with the database to be running.
    """
    driver_id = "test_driver_001"
    params = driver_db.get_driver_parameters(driver_id)

    assert params is not None
    assert params["driver_id"] == driver_id
    assert params["manufacturer"] == "Test Audio"
    assert params["model_name"] == "Model-X"
    assert params["fs_hz"] == 35.5
    assert params["qts"] == 0.35

def test_get_driver_parameters_not_found():
    """
    Test that retrieving a non-existent driver raises a ValueError.
    """
    with pytest.raises(ValueError, match="Driver with ID 'non_existent_driver' not found"):
        driver_db.get_driver_parameters("non_existent_driver") 