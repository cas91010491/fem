from thesis_sourcecode.src.model_training.patch_model_settings import PatchClassificationModel
from thesis_sourcecode.src.model_training.surface_points_model_settings import SurfacePointsModelForOnePatch
from thesis_sourcecode.src.model_training.signed_distance_model_settings import SignedDistanceModel
import numpy as np
from pdb import set_trace


if __name__ == "__main__":
    patch_classifier_name = "final_patch_model_edges-shape-512-512-bs-64"
    patch_classifier = PatchClassificationModel(name=patch_classifier_name)
    surface_points_regressor = []
    for i in range(96):
        name = f"final_surface_points_model_patch_{i}-shape-512-512-bs-64"
        print("sp name:" ,i)
        surface_points_regressor.append(SurfacePointsModelForOnePatch(name=name))
    signed_distance_regressor = SignedDistanceModel(name="testing_second_sd_model-shape-512-512-bs-64")

    points=[(0.1,0.2,0.4)]
    predicted_patches = patch_classifier.predict(points=points)
    # predicted_surface_points = []
    predicted_surface_points = np.zeros((len(points),2))
    for i in range(len(predicted_patches)):
        predicted_patch = predicted_patches[i]
        model_to_use = surface_points_regressor[predicted_patch]
        # predicted_surface_points.append(model_to_use.predict([points[i]]))
        predicted_surface_points[i]=(model_to_use.predict([points[i]]))
    predicted_signed_distance = signed_distance_regressor.predict(points=points)

    set_trace()

