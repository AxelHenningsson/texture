import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name="texture",
    version="0.0.0",
    author="Axel Henningsson",
    author_email="nilsaxelhenningsson@gmail.com",
    description="A collection of utility functions for analyzing crystallographic texture maps.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/AxelHenningsson/texture",
    project_urls={
        "Documentation": "https://axelhenningsson.github.io/texture/",
    },
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.8",
    install_requires=["matplotlib",
                      "numpy",
                      "scipy",
                      "pycifrw",
                      "xfab"]
)
