# Contributing to NetPharmPy

Thank you for your interest in contributing to NetPharmPy! We welcome contributions from the community to help make this network pharmacology analysis tool more robust and feature-rich.

## 🛠️ Development Setup

1.  **Clone the repository**
    ```bash
    git clone https://github.com/axp-knickei/netpharmpy.git
    cd netpharmpy
    ```

2.  **Create a virtual environment**
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate
    ```

3.  **Install dependencies**
    Install the package in editable mode along with required libraries:
    ```bash
    pip install -e .
    pip install -r requirements.txt
    pip install pytest pandas networkx pubchempy
    ```

## 🧪 Testing

We use `pytest` for testing. The test suite is divided into two categories:

### 1. Unit Tests (Fast & Offline)
Located in `tests/test_components.py`.
These tests mock external APIs (PubChem, STRING, Reactome) and file system operations. They are safe to run anywhere and verify the core logic of the data processing pipelines.

**Run unit tests:**
```bash
pytest tests/test_components.py
```

### 2. Integration Tests (Slower & Online)
Located in `tests/test_basic.py`.
These tests hit real external APIs. They ensure that our API wrappers are still compatible with the live services. Note that these may fail if the external services are down or rate-limited.

**Run all tests:**
```bash
pytest
```

## 📝 Coding Guidelines

*   **Style:** We follow standard Python PEP 8 conventions.
*   **Docstrings:** All public classes and methods must have docstrings (Google or NumPy style preferred) explaining arguments and return values.
*   **Type Hinting:** While not strictly enforced yet, type hints are encouraged for core logic.
*   **Imports:** Keep imports organized (Standard library > Third party > Local).

## 🚀 Pull Request Process

1.  Fork the repo and create your branch from `main`.
2.  If you've added code that should be tested, add tests.
3.  If you've changed APIs, update the documentation.
4.  Ensure the test suite passes.
5.  Make sure your code lints.
6.  Issue that pull request!

## 🐛 Reporting Bugs

If you find a bug, please open an issue on GitHub including:
*   The command you ran.
*   The expected output vs the actual output.
*   The log file content (found in your `outputs/` directory).
