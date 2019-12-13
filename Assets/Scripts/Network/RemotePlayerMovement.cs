﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RemotePlayerMovement : MonoBehaviour
{
    [SerializeField]
    private Vector3 nextPosition;
    private Quaternion nextRotation;
    
    public Vector3 NextPosition
    {
        set
        {
            nextPosition = value;
            StopCoroutine("moveAndRotate");
            StartCoroutine(moveAndRotate());
        }
    }

    public Quaternion NextRotation
    {
        set
        {
            nextRotation = value; 
            StopCoroutine("moveAndRotate");
            StartCoroutine(moveAndRotate());
        }
    }

    private float lerpInterval;
    // Start is called before the first frame update
    void Start()
    {
        lerpInterval = SimulationManager.GetInstance().MovementSendInterval;
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    IEnumerator moveAndRotate()
    {
        float elapsedTime = 0;
        Vector3 startPosition = transform.position;
        Quaternion startRotation = transform.rotation;
        while(elapsedTime < lerpInterval)
        {
            elapsedTime += Time.deltaTime;
            float delta = Mathf.Clamp(elapsedTime / lerpInterval, 0f, 0.99f);
            transform.position = Vector3.Lerp(startPosition, nextPosition, delta);
            if (startRotation != nextRotation)
            {
                transform.rotation = Quaternion.Lerp(startRotation, nextRotation, delta);
            }
            yield return new WaitForEndOfFrame();
        }
    }
}
