from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from alphafold import download_alphafold_pdb


class _FakeResponse:
    def __init__(self, *, json_data=None, content: bytes = b"") -> None:
        self._json = json_data
        self.content = content

    def raise_for_status(self) -> None:  # pragma: no cover - interface stub
        return None

    def json(self):  # pragma: no cover - interface stub
        return self._json


def test_download_alphafold_pdb_returns_absolute_path(tmp_path):
    output_file = tmp_path / "model.pdb"

    def _session_get(url, timeout=0):
        if "prediction" in url:
            return _FakeResponse(json_data=[{"pdbUrl": "https://example.org/model.pdb"}])
        return _FakeResponse(content=b"PDBDATA")

    with patch("alphafold.get_session") as mock_session:
        session = MagicMock()
        session.get.side_effect = _session_get
        mock_session.return_value = session

        result = download_alphafold_pdb("QTEST", str(output_file))

    assert result
    resolved = Path(result)
    assert resolved.is_file()
    assert resolved.read_bytes() == b"PDBDATA"
    assert resolved.is_absolute()


def test_download_alphafold_pdb_handles_missing_prediction(tmp_path):
    output_file = tmp_path / "model.pdb"

    def _session_get(url, timeout=0):
        return _FakeResponse(json_data=[])

    with patch("alphafold.get_session") as mock_session:
        session = MagicMock()
        session.get.side_effect = _session_get
        mock_session.return_value = session

        result = download_alphafold_pdb("QNONE", str(output_file))

    assert result == ""
    assert not output_file.exists()
