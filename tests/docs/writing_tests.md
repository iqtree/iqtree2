# Writing Tests

This document provides guidelines on how to write tests for the project using the provided test template.

## Using the Test Template

When writing a new test, you can use the provided `test_template.py` as a starting point. This template includes the necessary imports, a sample test function that utilizes the fixtures, and a `__main__` block that allows the test file to be run directly as a script during development of the test.

### Steps to Create a New Test Using the Template

1. **Copy the Template**: Copy the `test_template.py` file and rename it to something that reflects what you are testing. For example, if you are testing a feature called "example_feature", you might name your file `test_example_feature.py`.

2. **Update the Docstring**: Update the docstring at the beginning of the test function to reflect what your test is doing.

3. **Write the Test**: Replace the sample test logic in the `test_sample_feature` function with your own test logic to assert all the conditions that must be true for your test to pass.  keeping the following block so the test can be run directly.
```Python
if __name__ == "__main__":
    run_test_using(__file__)
```


4. **Run the Test**: Run the test directly using `python test_example_feature.py` to verify that it passes.

5. **Run All Tests**: Run all tests using `pytest` in the `\tests` root directory to verify that your new test does not break any existing tests.

6. **Auto-format the code**: Run `black .` to auto-format your code to ensure that it conforms to Python best practice style guidelines.

7. **Update requirements.txt**: If you have added any new Python dependencies to your test, update the `requirements.txt` file to include them.

7. **Commit and Push**: Commit your changes and push them to the repository.

### Fixtures
The following fixtures are available to all tests;
- **repo_paths**: returns a dictionary of paths to the project directories; 
  - repo_root: path to the repository root directory
  - build_dir: path to the build directory (location of iqtree binary)
  - tests_dir: path to the tests root directory
  - tests_data: path to any sample data that will be used in tests
- **temp_dir**: runs the test in a temporary directory. The temporary directory is automatically deleted after the test is run.  This is useful for testing features that create output files.
- **data_files**: takes a list of file names of files specified using a decorator `@pytest.mark.parametrize("data_files", [(["sample_file1.txt", "sample_file2.txt"])], indirect=True)` in the `/tests/data` directory and copies each file into the temporary directory that the test is run in.  This is useful for testing features that require input files.

### Parameterizing Tests

In some cases, you might want to run the same test function multiple times but with different sets of input data. The `pytest.mark.parametrize` decorator can be used for this purpose.

1. **Apply the Decorator**: Apply the `pytest.mark.parametrize` decorator above your test function. You need to specify the names of the arguments your test function will take, and a list of sets of data you want to pass to the test function.

    ```python
    @pytest.mark.parametrize("arg1, arg2", [(data1, data2), (data3, data4)])
    def test_example_feature(arg1, arg2):
        ...
    ```

3. **Write the Test Function**: Your test function should accept the arguments specified in the decorator and use them in the test logic.

    ```python
    def test_example_feature(arg1, arg2):
        # Test logic using arg1 and arg2
        ...
    ```