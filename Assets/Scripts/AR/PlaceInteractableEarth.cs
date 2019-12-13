using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.XR.ARFoundation;

[RequireComponent(typeof(ARTrackedImageManager))]
public class PlaceInteractableEarth : MonoBehaviour
{
    ARTrackedImageManager m_TrackedImageManager;
    public float scaleAdjust = 1.0f;
    private PlayerMovement localPlayerMovement;
    void Awake()
    {
        m_TrackedImageManager = GetComponent<ARTrackedImageManager>();
    }

    void OnEnable()
    {
        m_TrackedImageManager.trackedImagesChanged += OnTrackedImagesChanged;
    }

    void OnDisable()
    {
        m_TrackedImageManager.trackedImagesChanged -= OnTrackedImagesChanged;
    }
    void OnTrackedImagesChanged(ARTrackedImagesChangedEventArgs eventArgs)
    {
        foreach (var trackedImage in eventArgs.added)
        {
            trackedImage.transform.localScale = new Vector3(scaleAdjust, scaleAdjust, scaleAdjust);
            if (!localPlayerMovement)
            {
                localPlayerMovement = FindObjectOfType<PlayerMovement>();
            }

            if (localPlayerMovement)
            {
                localPlayerMovement.worldOrigin = trackedImage.transform;
            }
        }
    }
}
