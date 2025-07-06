from horn.data_ingestion import pdf_parser


def test_parse_pdf_datasheet():
    """
    Test that the parse_pdf_datasheet function exists and is callable.
    This is the first test in a TDD workflow.
    """
    assert callable(pdf_parser.parse_pdf_datasheet) 