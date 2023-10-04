def test_read_main(test_client) -> None:
    """Simple test of healthcheck entrypoint."""
    response = test_client.get("/healthcheck")
    assert response.status_code == 200
    assert response.json()["message"] == "API is running"
