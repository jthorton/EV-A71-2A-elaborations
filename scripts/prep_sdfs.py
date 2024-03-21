from rdkit import Chem
from pydantic import BaseModel, Field
from typing import Generator
from datetime import date
import pathlib


class Design(BaseModel):
    """A class to hold the selection information and metadata about this set of designs"""

    sdf_file: str = Field(..., description="The path to the sdf file")
    n_selected: int = Field(
        ...,
        description="The number of molecules which should be selected starting from the top"
        "we assume the sdf is sorted.",
    )
    ref_pdb: str = Field(
        ...,
        description="The fragalysis identifier of the pdb the ligands are posed against.",
    )
    ref_mols: list[str] = Field(
        ..., description="The list of reference ligands used to design this molecule."
    )
    rationale: str = Field(
        ..., description="The rationale used to design this molecule."
    )

    def molecules(self) -> Generator[Chem.Mol, None, None]:
        reader = Chem.SDMolSupplier(self.sdf_file, removeHs=False)
        return reader


class DesignDataset(BaseModel):
    """A class to hold information about a design dataset which can be used to prep the molecules for fragalysis"""

    ref_url: str = Field(
        ...,
        description="The url of the github readme with more information on the dataset.",
    )
    submitter_name: str = Field(
        ..., description="The name of the submitters of this dataset"
    )
    submitter_institution: str = Field(
        ..., description="The name of the institution of the submitter"
    )
    submitter_email: str = Field(..., description="The email of the submitter.")
    generation_date: str = Field(
        str(date.today()), description="The data the dataset was generated/ submitted."
    )
    method: str = Field(..., description="The method used to generate the dataset.")
    tags: list[str] = Field(
        ["cnnaffinityIC50", "cnnaffinity", "rationale"],
        description="The names of the extra tags which will be present in the molecules.",
    )
    designs: list[Design] = Field(
        [],
        description="The list of designs which should be combined to make our final selection.",
    )
    version: str = Field(
        "ver_1.2", description="The version of the fragalysis upload we are targeting."
    )


def workflow(dataset: DesignDataset, output_folder: str):
    """Prep the designdataset and create an SDF ready for fragalysis."""

    output_folder = pathlib.Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    sdf_file = output_folder.joinpath("fragalysis_ready.sdf")
    dataset_file = output_folder.joinpath("dataset_schema.json")
    with open(dataset_file.as_posix(), "w") as f:
        f.write(dataset.json(indent=2))

    with Chem.SDWriter(sdf_file.as_posix()) as writer:
        print("Building Dataset")
        # make the dummy molecule
        dummy_mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        dummy_mol.SetProp('_Name', dataset.version)
        data = dataset.dict(exclude={"version", "designs", "tags"})
        for k, v in data.items():
            dummy_mol.SetProp(k, v)
        # write the tags
        for tag in dataset.tags:
            dummy_mol.SetProp(tag, tag)
        # add the banner mol
        writer.write(dummy_mol)
        print("Dummy molecule added")

        # now for each dataset select the ligands and add the tags
        for design_set in dataset.designs:
            print(f"Preparing design set {design_set.sdf_file} selecting {design_set.n_selected}")
            selected = 0
            for molecule in design_set.molecules():
                # add the design set data
                molecule.SetProp("ref_pdb", design_set.ref_pdb)
                molecule.SetProp("ref_mols", ",".join(design_set.ref_mols))
                molecule.SetProp("rationale", design_set.rationale)
                name = molecule.GetProp("smiles").split(" ")[-1]
                molecule.SetProp("_Name", name)
                writer.write(molecule)
                selected += 1
                if selected == design_set.n_selected:
                    break

        print("Dataset built!")


if __name__ == "__main__":
    # build our dataset
    dataset = DesignDataset(
        ref_url="https://github.com/jthorton/EV-A71-2A-elaborations/tree/main/iteration2",
        submitter_name="JTHorton/DJCole",
        submitter_email="Josh.Horton@newcastle.ac.uk",
        submitter_institution="Newcastle University (UK)",
        method="FEGrow-Elaborations",
        designs=[
            Design(
                sdf_file="../iteration2/x0351_p2_elab/x0358_p2_available.sdf",
                n_selected=20,
                ref_pdb="Ax0310a",
                ref_mols=["Ax0351a"],
                rationale="Enamine real space elaborations around x0351."
            ),
            Design(
                sdf_file="../iteration2/x0365_elab/x0365_elab_available.sdf",
                n_selected=5,
                ref_pdb="Ax0310a",
                ref_mols=["Ax0365a"],
                rationale="Enamine real space elaborations around x0365."
            ),
            Design(
                sdf_file="../iteration2/x0528_x0556/x0528_elab_p2_p3_avilable.sdf",
                n_selected=20,
                ref_pdb="Ax0528a",
                ref_mols=["Ax0528a", "Ax0556a"],
                rationale="Enamine real space elaborations around a merge of x0528 and x0556."
            ),
            Design(
                sdf_file="../iteration2/x0875_elab/x0875_elab_available.sdf",
                n_selected=10,
                ref_pdb="Ax0875a",
                ref_mols=["Ax0875a"],
                rationale="Enamine real space elaborations around x0875."
            ),
            Design(
                sdf_file="../iteration2/x1019/x1019-elab.sdf",
                n_selected=5,
                ref_pdb="Ax0310a",
                ref_mols=["Ax1019a"],
                rationale="Enamine real space elaborations around x1019."
            ),
            Design(
                sdf_file="../iteration2/x1097_x0310/x1097-x0310_available.sdf",
                n_selected=20,
                ref_pdb="Ax0310a",
                ref_mols=["Ax1097a", "Ax0310a"],
                rationale="Enamine real space elaborations around a merge of x1097 and x0310."
            ),
            Design(
                sdf_file="../iteration2/x1097_x0739/x1097-x0739-quinoline-P3.sdf",
                n_selected=16,
                ref_pdb="Ax0310a",
                ref_mols=["Ax1097a", "Ax0739a"],
                rationale="Enamine real space elaborations around a merge of x1097 and x0739."
            ),
            Design(
                sdf_file="../iteration2/x1097_x0922/x1097-x0922_available.sdf",
                n_selected=10,
                ref_pdb="Ax0310a",
                ref_mols=["Ax1097a", "Ax0922a"],
                rationale="Enamine real space elaborations around a merge of x1097 and x0922."
            ),
        ]
    )

    workflow(dataset=dataset, output_folder="../iteration2/submission_1")
