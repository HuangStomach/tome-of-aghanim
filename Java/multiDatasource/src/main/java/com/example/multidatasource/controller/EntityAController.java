
package com.example.multidatasource.controller;

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

    @GetMapping
    public List<EntityA> getAll() {
        return repository.findAll();
    }

    @PostMapping
    public EntityA create(@RequestBody EntityA entity) {
        return repository.save(entity);
    }

    @GetMapping("/{id}")
    public EntityA getById(@PathVariable Long id) {
        return repository.findById(id).orElse(null);
    }

    @PutMapping("/{id}")
    public EntityA update(@PathVariable Long id, @RequestBody EntityA entity) {
        entity.setId(id);
        return repository.save(entity);
    }

    @DeleteMapping("/{id}")
    public void delete(@PathVariable Long id) {
        repository.deleteById(id);
    }
}
