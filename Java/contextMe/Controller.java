package com.example.contextMe.controller;

import com.example.multidatasource.model.a.EntityA;
import com.example.multidatasource.repository.a.EntityARepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@RestController
@RequestMapping("/api/a")
public class EntityAController {
    @Autowired
    private EntityARepository repository;
    @Autowired
    private Me me;

    @GetMapping()
    public ResponseEntity<?> list() {
        return ResponseEntity.ok().body(
            repository.findAllByUserId(me.getUserId()).getContent()
        );
    }
}
