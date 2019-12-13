using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PlayerMovement : MonoBehaviour
{
    private Vector3 lastPos;
    private Quaternion lastRot;
    private float lastSend = 0;
    NetworkController network;
    public Transform worldOrigin;
    private Transform defaultOrigin;
    private void Awake()
    {
        lastPos = transform.position;
        lastRot = transform.rotation;
        if (defaultOrigin == null)
        {
            defaultOrigin = Instantiate(new GameObject(), Vector3.zero, Quaternion.identity).transform;
        }

        if (worldOrigin == null)
        {
            Debug.Log("world origin unknown, using default");
            worldOrigin = Instantiate(new GameObject(), Vector3.zero, Quaternion.identity).transform;
        }
    }

    void FixedUpdate()
    {
        if (lastPos != transform.position || lastRot != transform.rotation)
        {
            // send update - no more frequently than once per second
            if (Time.time - SimulationManager.GetInstance().MovementSendInterval > lastSend)
            {
                // send movement!
                // Debug.Log("I'm sending!! " + transform.position);
                if (!network) network = FindObjectOfType<NetworkController>();
                Transform compareTransform = worldOrigin != null ? worldOrigin : defaultOrigin;
                network.HandleMovementUpdate(transform.position - compareTransform.position, transform.rotation, true);
                // update local comparators
                lastPos = transform.position;
                lastRot = transform.rotation;
                lastSend = Time.time;
            }
            
        }
    }
}
