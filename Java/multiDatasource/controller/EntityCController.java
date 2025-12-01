
package com.example.multidatasource.controller;

import com.example.multidatasource.model.c.EntityC;
import com.example.multidatasource.repository.c.EntityCRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@RestController
@RequestMapping("/api/c")
public class EntityCController {
    @Autowired
    private EntityCRepository repository;

    @GetMapping
    public List<EntityC> getAll() {
        return repository.findAll();
    }

    @GetMapping("/{id}")
    public EntityC getById(@PathVariable Long id) {
        return repository.findById(id).orElse(null);
    }
}
